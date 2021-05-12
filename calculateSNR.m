% (C) Copyright 2020 RnB FFT-analysis developers

function opt = calculateSNR(opt)
    % SNR analysis script adapted from
    % Xiaoqing Gao, Feb 27, 2020, Hangzhou xiaoqinggao@zju.edu.cn
    % dependent of CPP-BIDS, CPP-SPM and SPM functions

    %% set up experiment related info

    % number of steps per analysed period
    % Set default for RhythmFT and PitchFT analyses
    opt.whichHarmonics = [1, 2];

    % select which harmonics to take into account (always include 1st)
    % !!! careful: when we have slower frequency (4 steps per period),
    % don't select even harmonics in teh Block design. They overlap with
    % the sound-silence frequency (2 steps per period) !!!
    % Block step 4 use [1,3]
    % Block step 2 use [1,2]
    % FT step 4 use [1,2]
    if strcmpi(opt.taskName, 'RhythmBlock')
        if opt.nStepsPerPeriod == 2
            opt.whichHarmonics = [1, 2];
        elseif opt.nStepsPerPeriod == 4
            opt.whichHarmonics = [1, 3];
        end
    end

    % setup output directory
    for iSub = 1:numel(opt.subjects)

        subLabel = opt.subjects{iSub};
        
        destinationDir = createOutputDirectory(opt, subLabel);

        %% let's start

        % get mask image
        % use a predefined mask, only calculate voxels within the mask
        % below is same resolution as the functional images
        maskFileName = opt.funcMask{iSub};

        if opt.anatMask == 1
            maskFileName = opt.anatMaskFileName;
        end

        % load the mask
        maskHdr = spm_vol(maskFileName);
        maskImg = spm_read_vols(maskHdr);

        %     % output images
        %     outputHdr = spm_vol(outputImage);
        %     outputImg = spm_read_vols(outputHdr);

        % get functional image
        % we let SPM figure out what is in this BIDS data set
        opt = getSpecificBoldFiles(opt, subLabel);

        % add or count tot run number
        allRunFiles = opt.allFiles;

        %% setup parameters for FFT analysis
        % mri.repetition time(TR)
        repetitionTime = 1.75;

        % repetition of steps/categA
        patternDuration     = 12 * 0.190;
        segmentDuration     = 4 * patternDuration;
        stepDuration        = opt.nStepsPerPeriod * segmentDuration;

        % Number of vol before/after the rhythmic sequence (exp) are presented
        onsetDelay = 2;
        endDelay = 4;

        % use neighbouring 4 bins as noise frequencies
        cfg.binSize = 4;
        cfg.gap = 1;

        % set voxel and run numbers
        nVox = sum(maskImg(:) == 1);
        nRuns = length(allRunFiles);

        % number of samples (round to smallest even number)
        oldN = 105; % before resampling
        N = 104; % after resampling

        % calculate frequencies
        oddballFreq = 1 / stepDuration;
        oldFs =  1 / repetitionTime;
        fs = 1 / (182.4 / N);

        % frequencies
        freq = fs / 2 * linspace(0, 1, N / 2 + 1);

        % target frequency (this is a bad name, because it's a freq. bin INDEX)
        cfg.targetFrequency = round(N * oddballFreq / fs + 1);

        % harmonics of the target frequency and their bin indices
        cfg.harmonics = freq(cfg.targetFrequency) * opt.whichHarmonics;
        cfg.idxHarmonics = dsearchn(freq', cfg.harmonics');

        % number of bins for phase histogram
        cfg.histBin = 20;

        % threshold for choosing voxels for the phase distribution analysis
        cfg.thresh = 4;

        % preallocate to read all the runs and save
        allRunsBoldRaw = nan(N, nVox, nRuns);
        allRunsBoldDT = nan(N, nVox, nRuns);

        %%
        % design low-pass filter (to be 100% sure you prevent aliasing)
        fcutoff = fs / 4;
        transw  = .1;
        order   = round(7 * oldFs / fcutoff);
        shape   = [1 1 0 0];
        frex    = [0, fcutoff, fcutoff + fcutoff * transw, oldFs / 2] / (oldFs / 2);
        hz      = linspace(0, oldFs / 2, floor(oldN / 2) + 1);

        % get filter kernel
        filtkern = firls(order, frex, shape);

        % get kernel power spectrum
        filtkernX = abs(fft(filtkern, oldN)).^2;
        filtkernXdb = 10 * log10(abs(fft(filtkern, oldN)).^2);

        % % plot filter properties (visual check)
        %     figure
        %     plotedge = dsearchn(hz',fcutoff*3);
        %
        %     subplot(2,2,1)
        %     plot((-order/2:order/2)/fs,filtkern,'k','linew',3)
        %     xlabel('Time (s)')
        %     title('Filter kernel')
        %
        %     subplot(2,2,2), hold on
        %     plot(frex*fs/2,shape,'r','linew',1)
        %
        %     plot(hz,filtkernX(1:length(hz)),'k','linew',2)
        %     set(gca,'xlim',[0 fcutoff*3])
        %     xlabel('Frequency (Hz)'), ylabel('Gain')
        %     title('Filter kernel spectrum')
        %
        %     subplot(2,2,4)
        %     plot(hz,filtkernXdb(1:length(hz)),'k','linew',2)
        %     set(gca,'xlim',[0 fcutoff*3],'ylim',...
        %        [min([filtkernXdb(plotedge) filtkernXdb(plotedge)]) 5])
        %     xlabel('Frequency (Hz)'), ylabel('Gain')
        %     title('Filter kernel spectrum (dB)')

        %% Calculate SNR for each run
        % mypool = gcp('nocreate');
        % if isempty(mypool);
        %     parpool(4);
        % end

        %%
        % loop through runs
        for iRun = 1:nRuns

            fprintf('Read in file ... \n');

            % choose current BOLD file
            boldFile = allRunFiles{iRun};
            [~, boldFileName, ~] = fileparts(boldFile);

            % read/load bold file
            boldHdr = spm_vol(boldFile);
            boldImg = spm_read_vols(boldHdr);
            boldImg = reshape(boldImg, [size(boldImg, 1) * size(boldImg, 2) * ...
                                        size(boldImg, 3) size(boldImg, 4)]);

            % find cyclic volume
            totalVol = length(boldHdr);
            sequenceVol = totalVol - onsetDelay - endDelay;

            % remove the first 4 volumes, using this step to make the face stimulus onset at 0
            boldImg = boldImg(maskImg == 1, (onsetDelay + 1):(sequenceVol + onsetDelay));
            boldImg = boldImg';

            % control here, e.g. if the mask is bigger than the func img, we
            % would have zeros so we are getting rid off them below

            % save mask in the same size as boldImg and reallocate the matrix
            if iRun == 1

                % Getting rid off zeros
                zeroMask = all(boldImg == 0, 1);
                boldImg = boldImg(:, ~zeroMask);
                nVox = size(boldImg, 2);

                % remove also the mask' image zeros
                removeMaskZeros(maskFileName, boldFile);

                % reload the mask
                maskHdr = spm_vol(maskFileName);
                maskImg = spm_read_vols(maskHdr);

                % reload the new mask

                % reallocate with the correct voxel size
                allRunsBoldRaw = nan(N, nVox, nRuns);
                allRunsBoldDT = nan(N, nVox, nRuns);
            end

            % filter and interpolate
            boldResampled = zeros(N, size(boldImg, 2));

            % interpolate (resample)
            timeVector = linspace(1, oldN, N);
            oldTimeVector = 1:oldN;

            for iVoxel = 1:size(boldImg, 2)
                % low-pass filter
                boldFilter = filtfilt(filtkern, 1, boldImg(:, iVoxel));
                % interpolate
                boldResampled(:, iVoxel) = interp1(oldTimeVector, boldFilter, ...
                                                   timeVector, 'spline');
            end

            % check that sizes ale okay
            if size(boldImg, 1) ~= oldN
                error('check number of volumes for the original');
            end
            if size(boldResampled, 1) ~= N
                error('check number of volumes for the original');
            end

            % remove linear trend
            boldDetrend = detrend(boldResampled);

            allRunsBoldRaw(:, :, iRun) = boldResampled;
            allRunsBoldDT(:, :, iRun) = boldDetrend;

            % save the output
            if opt.saveEachRun == 1

                fprintf('Saving each run output... \n');
                
                % run FFT
                [targetZ, cfg, ~] = calculateFourier(boldDetrend, boldResampled, cfg);

                %     % unused parameters for now
                %     targetPhase = cfg.targetPhase;
                %     targetSNRsigned = cfg.targetSNRsigned;
                %     tSNR = cfg.tSNR;
                %     %

                newFileName = ['SNR_', boldFileName, '.nii'];
                writeMap(targetZ, maskHdr, maskImg, newFileName, destinationDir);

            end
        end

        %% Calculate SNR for the averaged time course of the all runs
        fprintf('Calculating average... \n');

        % rename the file name for saving
        boldFileName = regexprep(boldFileName, 'run-(\d*)_', '');

        % average runs in the time domain
        avgBold = mean(allRunsBoldDT, 3);
        avgRawBold = mean(allRunsBoldRaw, 3);

        %% run FFT on average bold
        % ----------------------------------------------------------------
        % Zscore at target frequency. 
        % (1) calculate magnitude specta using FFT
        % (2) exctract magnitude at the target frequency bin -> this is
        % "signal"
        % (3) extract magnitudes at surrounding frequency bins -> this is noise
        % (4) calculate zscore as (signal-mean(noise)) / std(noise)
        % do this for each voxel and return in array AvgZTarget
        [AvgZTarget, cfg, FT] = calculateFourier(avgBold, avgRawBold, cfg);
        
        % save map as nii
        newFileName = ['AvgZTarget_', boldFileName, '.nii'];
        zmapImg = writeMap(AvgZTarget, maskHdr, maskImg, newFileName, destinationDir);

        %% subtraction SNR 
        % baseline-correct magnitude spectra (subtracting surrounding bins)
        
        % convert from complex FFT output to magnitudes 
        mX = abs(FT); 
        
        % define the baseline correction function (we will simply subtract the
        % mean magnitude at surrounding bins from the signal, i.e. center bin)
        blfun = @(signal, noise) signal - noise;
        % go over frequency bins, and apply the function to each bin in the
        % spectra
        mXbl = baselineCorrect(mX, cfg, 'fun', blfun);
        
        % plot best voxels
        voxelNbToPlot = 10;
        coordTargetFreq = getVoxelCoordinate(boldHdr, zmapImg, maskImg, voxelNbToPlot);
        % for a clear representation, does it make sense to use fixed scale
        % across voxels?
        f = plotmXBestVox(freq, mXbl, coordTargetFreq, cfg.idxHarmonics);

        % save figure
        newFileName = 'AvgZTarget-bestVox_';
        fileName = fullfile(destinationDir, [newFileName, boldFileName, '.fig']);
        saveas(f, fileName);
        close(f);

        % ----------------------------------------------------------------
        % Zscore at target frequency and harmonics. We will first pool
        % (average) magnitudes across the different harmonics. To this end, we
        % will cut segments around the target frequency and each harmonic from
        % the long spectra. We will average these segments and from this
        % avergae, we will extract the zscore (signal vs. noise). 
        
        % cut out segments of frequency bins around each frequency of interest
        % (target and its harmonics) and average these segments 
        mXavgHarmonics = getAvgHarmonics(mX, cfg.idxHarmonics, cfg.binSize);

        % define indices of positions in the averaged segment that represent
        % noise (note that the center position is the frequency bin where we
        % expect signal, and noise is on either side of this bin) 
        noiseFs = [(cfg.binSize + 1 - cfg.binSize / 2 - cfg.gap) ...
                   :(cfg.binSize + 1 - 1 - cfg.gap) ...
                   (cfg.binSize + 1 + 1 + cfg.gap) ...
                   :(cfg.binSize + 1 + cfg.binSize / 2 + cfg.gap)];

        % extract magnitudes of signal and noise from the average segment 
        AmpNoise = mXavgHarmonics(noiseFs, :);
        NoiseMean = mean(AmpNoise, 1);
        NoiseSD = std(AmpNoise, 0, 1);
        % get zscore 
        targetHarmonicsZ = (mXavgHarmonics(cfg.binSize + 1, :) - NoiseMean) ./ NoiseSD;

        % save map as nii
        newFileName = ['AvgZHarmonics_', boldFileName, '.nii'];
        zHarmonicsImg = writeMap(targetHarmonicsZ, maskHdr, maskImg, newFileName, destinationDir);

        coordHarmonics = getVoxelCoordinate(boldHdr, zHarmonicsImg, maskImg, voxelNbToPlot);
        
        % plot averge segment from the spectra for best voxels
        % for a clear representation, does it make sense to use fixed scale
        % across voxels?
        f = plotAvgHarmBestVox(mXavgHarmonics, coordHarmonics); 
        
        % save figure
        newFileName = ['AvgZHarmonics-bestVoxMean_', boldFileName, '.fig'];
        saveas(f, fullfile(destinationDir, newFileName));
        close(f);

        % plot whole spectra for best voxels and save figure
        % the coordinate of best harmonics differ from coord of highest SNR
        % on target freuency - for Fig1 differs from Fig3
        f = plotmXBestVox(freq, mXbl, coordHarmonics, cfg.idxHarmonics);
        newFileName = ['AvgZHarmonics-bestVox_', boldFileName, '.fig'];
        saveas(f, fullfile(destinationDir, newFileName));
        close(f);

        %% ratio SNR 
        % baseline-correct magnitude spectra (taking ration with magnitudes
        % at surrounding bins)
        
        % define baseline correction function (we will simply take a ratio
        % between the mean magnitude the signal frequency bin, and mean magnitude
        % at surrounding bins (from both sides) 
        blfun = @(x, y) x ./ y;
        % go over frequency bins, and apply the function to each bin in the
        % spectra
        mXblRatio = baselineCorrect(abs(FT), cfg, 'fun', blfun);
        
        % extract the value at target frequency (separately for each voxel) 
        targetRatio = mXblRatio(cfg.targetFrequency, :);

        % save map as nii
        newFileName = ['AvgRatioTarget_', boldFileName, '.nii'];
        targetRatioImg = writeMap(targetRatio, maskHdr, maskImg, newFileName, destinationDir);

        % plot best voxels
        coordRatio = getVoxelCoordinate(boldHdr, targetRatioImg, maskImg, voxelNbToPlot);        
        
        f = plotmXBestVox(freq, mXblRatio, coordRatio, cfg.idxHarmonics, 'ratio');
        
        % save figure
        newFileName = ['AvgRatioTarget-bestVox_', boldFileName, '.fig'];
        saveas(f, fullfile(destinationDir, newFileName));
        close(f);

    end
end

function zmapImg = writeMap(data2write, maskHdr, maskImg, newFileName, destinationDir)

    % create template hdr to be saved
    zmapHdr = maskHdr;
    zmapHdr.fname = spm_file(zmapHdr.fname, 'path', destinationDir);
    zmapHdr.fname = spm_file(zmapHdr.fname, 'filename', newFileName);

    % get dimensions & allocate 3-D img
    zmapImg = zeros(zmapHdr.dim);

    % get mask index for non-zero values &
    % assign z-scores from 1-D to their correcponding 3-D location
    zmapImg(find(maskImg > 0)) = data2write; %#ok<FNDSB>

    % save result as .nii file
    spm_write_vol(zmapHdr, zmapImg);

end

function removeMaskZeros(maskFileName, boldFile)

    % read mask image
    maskHdr = spm_vol(maskFileName);
    maskImg = spm_read_vols(maskHdr);

    % read bold image
    boldHdr = spm_vol(boldFile);
    boldImg = spm_read_vols(boldHdr);

    % take 1 image from whole run
    boldImg1 = boldImg(:, :, :, 1);

    % make bold image binary
    boldImg1(boldImg1 ~= 0) = 1;

    % overlap (temp = bold + mask )
    temp = maskImg + boldImg1;
    temp(temp == 1) = 0;
    temp(temp == 2) = 1;

    % save new mask
    spm_write_vol(maskHdr, temp);

end

% function writeMap(data2write, maskFileName, newFileName, destinationDir)
%
%
%     % write map of extracted values
%     mask_new = load_untouch_nii(maskFileName);
%
%     zmap3Dmask = zeros(size(mask_new.img));
%     zmap3Dmask(find(mask_new.img > 0)) = data2write;
%
%     new_nii = make_nii(zmap3Dmask);
%     new_nii.hdr = mask_new.hdr;
%
%     dims = size(mask_new.img);
%     new_nii.hdr.dime.dim(2:5) = [dims(1) dims(2) dims(3) 1];
%
%
%     save_nii(new_nii, fullfile(destinationDir,newFileName));
%
% end
