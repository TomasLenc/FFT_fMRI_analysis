% (C) Copyright 2020 RnB FFT-analysis developers

function destinationDir = createOutputDirectory(opt)
% sets a destination directory for FFT analysis results

  subjectDestDir = fullfile(opt.derivativesDir, '..', 'FFT_RnB_funcmask');
  if opt.anatMask
    subjectDestDir = fullfile(opt.derivativesDir, '..', 'FFT_RnB_anatmask');
  end

  subject = ['sub-', opt.subjects{1}];
  session = 'ses-001';
  stepFolder = ['step', num2str(opt.stepSize)];
  dirsToMake = {subject, session, stepFolder};

  % create subject folder witn subfolders if doesn't exist
  if ~exist(fullfile(subjectDestDir, subject, session, stepFolder), 'dir')
    for idir = 1:length(dirsToMake)
      Thisdir = fullfile(subjectDestDir, dirsToMake{1:idir});
      if ~exist(Thisdir)
        mkdir(Thisdir);
      end
    end
  end

  % output the results
  destinationDir =  fullfile(subjectDestDir, subject, session, stepFolder);

end