% (C) Copyright 2020 RnB FFT-analysis developers

function opt = getSpecificBoldFiles(opt)
% gets the specificed bold files for the FFT analysis
% dependend on cpp-spm and bids-matlab functions. 


  % we let SPM figure out what is in this BIDS data set
  [~, opt, BIDS] = getData(opt);

  subID = opt.subjects(1);

  %% Get functional files for FFT
  % identify sessions for this subject
  [sessions, nbSessions] = getInfo(BIDS, subID, opt, 'Sessions');

  % get prefix for smoothed image
  [prefix, ~] = getPrefix('ffx', opt, opt.FWHM);

  allFiles = [];
  sesCounter = 1;

  for iSes = 1:nbSessions        % For each session

    % get all runs for that subject across all sessions
    [runs, nbRuns] = getInfo(BIDS, subID, opt, 'Runs', sessions{iSes});

    for iRun = 1:nbRuns

      % get the filename for this bold run for this task
      [fileName, subFuncDataDir] = getBoldFilename( ...
                                                   BIDS, ...
                                                   subID, sessions{iSes}, ...
                                                   runs{iRun}, opt);

      % check that the file with the right prefix exist
      files = validationInputFile(subFuncDataDir, fileName, prefix);

      % add the files to list
      allFilesTemp = cellstr(files);
      allFiles = [allFiles; allFilesTemp]; %#ok<AGROW>
      sesCounter = sesCounter + 1;

    end
  end

  opt.allFiles = allFiles;

  %% get the masks for FFT

  % get mean image
  [meanImage, meanFuncDir] = getMeanFuncFilename(BIDS, subID, opt);
  meanFuncFileName = fullfile(meanFuncDir, meanImage);

  % normalized image option by adding prefix w-
  if strcmp(opt.space, 'MNI')
    meanFuncFileName = fullfile(meanFuncDir, ['w', meanImage]);
  end

  % think about it again % % % %
  % instead of segmented meanfunc image here
  % get native-spaced resliced anat (cpp-spm pipeline) image:
  [~, meanImageName, ext] = fileparts(meanImage);
  anatMaskFileName = fullfile(meanFuncDir, ...
                              [meanImageName, '_mask', ext]);

  opt.anatMaskFileName = anatMaskFileName;
  opt.funcMaskFileName = meanFuncFileName;

  % save prefix
  opt.prefix = prefix;

end