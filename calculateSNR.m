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

    %% let's start
    % we let SPM figure out what is in this BIDS data set
    opt = getSpecificBoldFiles(opt);

    % add or count tot run number
    allRunFiles = opt.allFiles;

    % use a predefined mask, only calculate voxels within the mask
    maskFileName = opt.anatMaskFileName;

    if ~opt.anatMask == 1
        maskFileName = makeFuncIndivMask(opt);
    end

    maskFile = spm_vol(maskFileName);
    mask = spm_read_vols(maskFile);

    %% setup parameters for FFT analysis
    % mri.repetition time(TR)
    repetitionTime = 1.75;

    % repetition of steps/categA
    patternDuration     = 12 * 0.190;
    segmentDuration     = 4 * patternDuration;
    stepDuration        = opt.nStepsPerPeriod * segmentDuration;

    % setup output directory
    opt.destinationDir = createOutputDirectory(opt);

    % Number of vol before/after the rhythmic sequence (exp) are presented
    onsetDelay = 2;
    endDelay = 4;

    % use neighbouring 4 bins as noise frequencies
    cfg.binSize = 4;
    cfg.gap = 1;

    % set voxel and run numbers
    RunPattern = struct();
    nVox = sum(mask(:) == 1);
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
        boldFileName = allRunFiles{iRun};

        % read/load bold file
        boldFile = spm_vol(boldFileName);
        signal = spm_read_vols(boldFile); % check the load_untouch_nii to compare
        signal = reshape(signal, [size(signal, 1) * size(signal, 2) * ...
                                  size(signal, 3) size(signal, 4)]);

        % find cyclic volume
        totalVol = length(spm_vol(boldFileName));
        sequenceVol = totalVol - onsetDelay - endDelay;

        % remove the first 4 volumes, using this step to make the face stimulus onset at 0
        Pattern = signal(mask == 1, (onsetDelay + 1):(sequenceVol + onsetDelay));
        Pattern = Pattern';

        % interpolate (resample)
        xi = linspace(0, oldN, N);

        % filter and interpolate
        patternResampled = zeros(N, size(Pattern, 2));

        for voxi = 1:size(Pattern, 2)
            % low-pass filter
            PatternFilt = filtfilt(filtkern, 1, Pattern(:, voxi));
            % interpolate
            patternResampled(:, voxi) = interp1([1:oldN], PatternFilt, xi, 'spline');
        end

        % check that sizes ale okay
        if size(Pattern, 1) ~= oldN
            error('check number of volumes for the original');
        end
        if size(patternResampled, 1) ~= N
            error('check number of volumes for the original');
        end

        % remove linear trend
        patternDetrend = detrend(patternResampled);

        [targetZ, cfg] = calculateFourier(patternDetrend, patternResampled, cfg);

        %     % unused parameters for now
        %     targetPhase = cfg.targetPhase;
        %     targetSNRsigned = cfg.targetSNRsigned;
        %     tSNR = cfg.tSNR;
        %     %

        allRunsRaw(:, :, iRun) = patternResampled;
        allRunsDT(:, :, iRun) = patternDetrend;

        fprintf('Saving ... \n');

        % z-scored 1-D vector
        zmapmasked = targetZ;

        % allocate 3-D img
        % get the mask
        mask_new = load_untouch_nii(maskFileName);
        zmap3Dmask = zeros(size(mask_new.img));

        % get mask index
        maskIndex = find(mask_new.img == 1);

        % assign z-scores from 1-D to their correcponding 3-D location
        zmap3Dmask(maskIndex) = zmapmasked;

        new_nii = make_nii(zmap3Dmask);

        new_nii.hdr = mask_new.hdr;

        % get dimensions to save
        dims = size(mask_new.img);
        new_nii.hdr.dime.dim(2:5) = [dims(1) dims(2) dims(3) 1];

        %   % save the results
        %   FileName = fullfile(opt.destinationDir, ['SNR_', boldFileName, ext]);
        %   save_nii(new_nii, FileName);

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
    fileName = [getOutputFileName(valueName, boldFileName, opt), '.fig'];
    saveas(f, fileName);
    close(f);

    % save map as nii
    valueName = 'AvgZTarget';
    writeMap(targetZ, maskFileName, valueName, opt);

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
    fileName = [getOutputFileName(valueName, boldFileName, opt), '.fig'];
    saveas(f, fileName);
    close(f);

    % plot best voxels and save figure
    f = plotmXBestVox(freq, mXSNRAmp, targetHarmonicsZ, 10, cfg.idxHarmonics);
    valueName = 'AvgZHarmonics-bestVox';
    fileName = [getOutputFileName(valueName, boldFileName, opt), '.fig'];
    saveas(f, fileName);
    close(f);

    % save map as nii
    valueName = 'AvgZHarmonics';
    writeMap(targetHarmonicsZ, maskFileName, valueName, opt);

    % ----------------------------------------------------------------
    % SNR at target frequency as ratio to neighbouring bins
    blfun = @(x, y) x ./ y;
    mXSNRRatio = baselineCorrect(abs(FT), cfg, 'fun', blfun);
    targetSNRRatio = mXSNRRatio(cfg.targetFrequency, :);

    % plot best voxels
    f = plotmXBestVox(freq, mXSNRRatio, targetSNRRatio, 10, cfg.idxHarmonics, 'ratio');

    % save figure
    valueName = 'AvgRatioTarget-bestVox';
    fileName = [getOutputFileName(valueName, boldFileName, opt), '.fig'];
    saveas(f, fileName);
    close(f);

    % save map as nii
    valueName = 'AvgRatioTarget';
    writeMap(targetSNRRatio, maskFileName, valueName, opt);

end

function writeMap(data2write, maskFileName, valueName, opt)

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
    fileName = [getOutputFileName(valueName, boldFileName, opt), '.nii'];

    % save output file
    save_nii(new_nii, fileName);

end

function fileName = getOutputFileName(valueName, boldFileName, opt)

    fileName = fullfile(opt.destinationDir, [valueName, boldFileName]);
end

% % set save filename and save results as .nii
% FileName = fullfile(opt.destinationDir, ['AvgSNR_', boldFileName, ext]);
% save_nii(new_nii, FileName);
