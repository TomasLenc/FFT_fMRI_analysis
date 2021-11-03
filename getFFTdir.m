% (C) Copyright 2020 RnB FFT-analysis developers

function fftDir = getFFTdir(opt, subLabel)
    % sets a destination directory for FFT analysis results

    maskType = 'func';

    mainDir = fullfile(opt.derivativesDir, '..', 'rnb_fft');

    % in the future omit hardcoding of subject - loop through instead
    subject = ['sub-', subLabel];
    session = 'ses-001';
    fftDirName = ['space-', opt.space, ...
                     '_FWHM-', num2str(opt.FWHM), ...
                     '_mask-', maskType, ...
                     '_step-', num2str(opt.nStepsPerPeriod)];

    % output the results
    fftDir =  fullfile(mainDir, subject, session, fftDirName);

    % create subject folder witn subfolders if doesn't exist

    if ~exist(fftDir, 'dir')

        dirsToMake = {subject, session, fftDirName};

        for idir = 1:length(dirsToMake)

            Thisdir = fullfile(mainDir, dirsToMake{1:idir});

            if ~exist(Thisdir)
                mkdir(Thisdir);
            end

        end

    end

end
