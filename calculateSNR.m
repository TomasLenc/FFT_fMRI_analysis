% (C) Copyright 2020 RnB FFT-analysis developers

function opt = calculateSNR(opt)
    % SNR analysis script adapted from
    % Xiaoqing Gao, Feb 27, 2020, Hangzhou xiaoqinggao@zju.edu.cn
    % dependent of CPP-BIDS, CPP-SPM and SPM functions

    %% set up experiment related info

    % number of steps per analysed period
    opt.nStepsPerPeriod = 4;

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
    else
        strcmpi(opt.taskName, 'RhythmFT');
        if opt.nStepsPerPeriod == 4
            opt.whichHarmonics = [1, 2];
        end
    end

    % setup output directory
    destinationDir = createOutputDirectory(opt);
    
    %% let's start
    
    % get mask image
    % use a predefined mask, only calculate voxels within the mask
    % below is same resolution as the functional images
    maskFileName = opt.funcMask;
    
    if opt.anatMask == 1
        maskFileName = opt.anatMaskFileName;
    end

    % load the mask
    maskHdr = spm_vol(maskFileName);
    maskImg = spm_read_vols(maskHdr);

    % get functional image
    % we let SPM figure out what is in this BIDS data set
    opt = getSpecificBoldFiles(opt);

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
    %RunPattern = struct();
    nVox = sum(maskImg(:) > 0);
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
    allRunsRaw = nan(N, nVox, nRuns);
    allRunsDT = nan(N, nVox, nRuns);

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
    % figure
    % plotedge = dsearchn(hz',fcutoff*3);
    %
    % subplot(2,2,1)
    % plot((-order/2:order/2)/fs,filtkern,'k','linew',3)
    % xlabel('Time (s)')
    % title('Filter kernel')
    %
    % subplot(2,2,2), hold on
    % plot(frex*fs/2,shape,'r','linew',1)
    %
    % plot(hz,filtkernX(1:length(hz)),'k','linew',2)
    % set(gca,'xlim',[0 fcutoff*3])
    % xlabel('Frequency (Hz)'), ylabel('Gain')
    % title('Filter kernel spectrum')
    %
    % subplot(2,2,4)
    % plot(hz,filtkernXdb(1:length(hz)),'k','linew',2)
    % set(gca,'xlim',[0 fcutoff*3],'ylim',...
    %    [min([filtkernXdb(plotedge) filtkernXdb(plotedge)]) 5])
    % xlabel('Frequency (Hz)'), ylabel('Gain')
    % title('Filter kernel spectrum (dB)')

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
        [~,boldFileName, ext] = fileparts(boldFile);
        
        % read/load bold file
        boldHdr = spm_vol(boldFile);
        boldImg = spm_read_vols(boldHdr); 
        boldImg = reshape(boldImg, [size(boldImg, 1) * size(boldImg, 2) * ...
                                  size(boldImg, 3) size(boldImg, 4)]);

        % find cyclic volume
        totalVol = length(boldHdr);
        sequenceVol = totalVol - onsetDelay - endDelay;

        % remove the first 4 volumes, using this step to make the face stimulus onset at 0
        boldImg = boldImg(maskImg > 0, (onsetDelay + 1):(sequenceVol + onsetDelay));
        boldImg = boldImg';
        
        % filter and interpolate
        boldResampled = zeros(N, size(boldImg, 2));
        
        % interpolate (resample)
        timeVector = linspace(1, oldN, N);
        oldTimeVector = 1:oldN;
        

        for iVoxel = 1:size(boldImg, 2)
            % low-pass filter
            boldFilter = filtfilt(filtkern, 1, boldImg(:, iVoxel));
            % interpolate
            boldResampled(:, iVoxel) = interp1(oldTimeVector, boldFilter,...
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

        % run FFT
        [targetZ, cfg] = calculateFourier(boldDetrend, boldResampled, cfg);

        %     % unused parameters for now
        %     targetPhase = cfg.targetPhase;
        %     targetSNRsigned = cfg.targetSNRsigned;
        %     tSNR = cfg.tSNR;
        %     %

        allRunsRaw(:, :, iRun) = boldResampled;
        allRunsDT(:, :, iRun) = boldDetrend;

        fprintf('Saving ... \n');

        % create template hdr to be saved
        zmapHdr = maskHdr;
        zmapnewName = fullfile(destinationDir, ...
                               ['SNR_', boldFileName, ext]);
        %zmapHdr.fname =  zmapnewName;
        zmapHdr.fname = spm_file(zmapHdr.fname,'filename',zmapnewName);
        
        % % % 
        
        % use spm_file with another option than "filename" to change the
        % saving path. at the moment it still saves into :
        % '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/PitchFT/derivatives/cpp_spm/SNR_s3wuasub-011_ses-001_task-PitchFT_run-001_bold.nii'
        % fpath did not work
        
        % % %
        % get dimensions & allocate 3-D img
        zmapImg = zeros(zmapHdr.dim);

        % get mask index for non-zero values &
        % assign z-scores from 1-D to their correcponding 3-D location
        zmapImg(find(maskImg > 0)) = targetZ; %#ok<FNDSB>

        % save result as .nii file
        spm_write_vol(zmapHdr, zmapImg);

    end

    %% Calculate SNR for the averaged time course of the all runs
    fprintf('Calculating average... \n');

    % average runs in the time domain
    avgPattern = mean(allRunsDT, 3);
    avgrawPattern = mean(allRunsRaw, 3);

    % ----------------------------------------------------------------
    % zscore at target frequency

    blfun = @(x, y) x - y;
    [targetZ, cfg, FT] = calculateFourier(avgPattern, avgrawPattern, cfg);
    mXSNRAmp = baselineCorrect(abs(FT), cfg, 'fun', blfun);

    % plot best voxels
    f = plotmXBestVox(freq, mXSNRAmp, targetZ, 10, cfg.idxHarmonics);

    % save figure
    valueName = 'AvgZTarget-bestVox';
    fileName = [getOutputFileName(valueName, boldFile, destinationDir), '.fig'];
    saveas(f, fileName);
    close(f);

    % save map as nii
    valueName = 'AvgZTarget';
    writeMap(targetZ, maskFileName, valueName, destinationDir);

    % ----------------------------------------------------------------
    % zscore at target frequency and harmonics (average amp first)
    mXavgHarmonics = getAvgHarmonics(abs(FT), cfg.idxHarmonics, cfg.binSize);

    noiseFs = [(cfg.binSize + 1 - cfg.binSize / 2 - cfg.gap) ...
               :(cfg.binSize + 1 - 1 - cfg.gap) ...
               (cfg.binSize + 1 + 1 + cfg.gap) ...
               :(cfg.binSize + 1 + cfg.binSize / 2 + cfg.gap)];

    AmpNoise = mXavgHarmonics(noiseFs, :);
    NoiseMean = mean(AmpNoise, 1);
    NoiseSD = std(AmpNoise, 0, 1);
    targetHarmonicsZ = (mXavgHarmonics(cfg.binSize + 1, :) - NoiseMean) ./ NoiseSD;

    % plot best voxels
    [~, idxSorted] = sort(targetHarmonicsZ, 'descend');
    idxSorted(isnan(targetHarmonicsZ(idxSorted))) = [];
    nMax = 10;
    f = figure('color', 'white', 'Position', [131 728 1744 140]);
    pnl = panel(f);
    pnl.pack('h', nMax);
    for iVox = 1:nMax
        mXbest = mXavgHarmonics(:, idxSorted(iVox));
        pnl(iVox).select();
        h = stem(mXbest, ...
                 'marker', 'none', ...
                 'color', [0.6, 0.6, 0.6], ...
                 'linew', 4);
        hold on;
        h = stem(cfg.binSize + 1, mXbest(cfg.binSize + 1), ...
                 'marker', 'none', ...
                 'color', 'r', ...
                 'linew', 4);
        title(sprintf('z=%.2f  vox=%d', targetHarmonicsZ(idxSorted(iVox)), ...
                      idxSorted(iVox)));
        set(gca, 'xtick', []);
    end
    pnl.ylabel('amplitude');
    pnl.marginbottom = 3;

    % save figure
    valueName = 'AvgZHarmonics-bestVoxMean';
    fileName = [getOutputFileName(valueName, boldFile, destinationDir), '.fig'];
    saveas(f, fileName);
    close(f);

    % plot best voxels and save figure
    f = plotmXBestVox(freq, mXSNRAmp, targetHarmonicsZ, 10, cfg.idxHarmonics);
    valueName = 'AvgZHarmonics-bestVox';
    fileName = [getOutputFileName(valueName, boldFile, destinationDir), '.fig'];
    saveas(f, fileName);
    close(f);

    % save map as nii
    valueName = 'AvgZHarmonics';
    writeMap(targetHarmonicsZ, maskFileName, valueName, destinationDir);

    % ----------------------------------------------------------------
    % SNR at target frequency as ratio to neighbouring bins
    blfun = @(x, y) x ./ y;
    mXSNRRatio = baselineCorrect(abs(FT), cfg, 'fun', blfun);
    targetSNRRatio = mXSNRRatio(cfg.targetFrequency, :);

    % plot best voxels
    f = plotmXBestVox(freq, mXSNRRatio, targetSNRRatio, 10, cfg.idxHarmonics, 'ratio');

    % save figure
    valueName = 'AvgRatioTarget-bestVox';
    fileName = [getOutputFileName(valueName, boldFile, destinationDir), '.fig'];
    saveas(f, fileName);
    close(f);

    % save map as nii
    valueName = 'AvgRatioTarget';
    writeMap(targetSNRRatio, maskFileName, valueName, destinationDir);

end

function writeMap(data2write, maskFileName, valueName, destinationDir)

    % write map of extracted values
    mask_new = load_untouch_nii(maskFileName);

    maskIndex = find(mask_new.img == 1);

    dims = size(mask_new.img);

    zmapmasked = data2write;

    zmap3Dmask = zeros(size(mask_new.img));

    zmap3Dmask(maskIndex) = zmapmasked;

    new_nii = make_nii(zmap3Dmask);

    new_nii.hdr = mask_new.hdr;

    new_nii.hdr.dime.dim(2:5) = [dims(1) dims(2) dims(3) 1];

    % create a filename
    fileName = [getOutputFileName(valueName, boldFileName, destinationDir), '.nii'];

    % save output file
    save_nii(new_nii, fileName);

end

function fileName = getOutputFileName(valueName, boldFileName, destinationDir)

    fileName = fullfile(destinationDir, [valueName, boldFileName]);
end

% % set save filename and save results as .nii
% FileName = fullfile(opt.destinationDir, ['AvgSNR_', boldFileName, ext]);
% save_nii(new_nii, FileName);
