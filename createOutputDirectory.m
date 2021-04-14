% (C) Copyright 2020 RnB FFT-analysis developers

function destinationDir = createOutputDirectory(opt, subLabel)
    % sets a destination directory for FFT analysis results

    maskType = 'func';
    if opt.anatMask 
        maskType = 'anat';
    end
    
    fftDir = fullfile(opt.derivativesDir, '..', 'rnb_fft');


    
    % in the future omit hardcoding of subject - loop through instead
    subject = ['sub-', subLabel];
    session = 'ses-001';
    subjectFftDir = ['space-', opt.space, ...
                             '_FWHM-', num2str(opt.FWHM), ...
                             '_mask-', maskType, ...
                             '_step-', num2str(opt.nStepsPerPeriod)];

    

    % output the results
    destinationDir =  fullfile(fftDir, subject, session, subjectFftDir);
    
    % create subject folder witn subfolders if doesn't exist
    
    if ~exist(destinationDir, 'dir')
        
        dirsToMake = {subject, session, subjectFftDir};
        
        for idir = 1:length(dirsToMake)
            
            Thisdir = fullfile(fftDir, dirsToMake{1:idir});
            
            if ~exist(Thisdir)
                mkdir(Thisdir);
            end
            
        end
        
        
    end

    

end
