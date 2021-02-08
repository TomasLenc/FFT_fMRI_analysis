% (C) Copyright 2020 RnB FFT-analysis developers

function mask = makeFuncIndivMask(opt)

  % function cpp-spm repo/functions and create skull stripped 
  % mean functional image

  % read the mask image
  image = opt.funcMaskFileName;
  [imagePath, imageName, ext] = fileparts(image);


  %SPM skull stripping - with Anat atm
  [~, opt, BIDS] = getData(opt);
    subID = opt.subjects{1};
  
  opt.orderBatches.segment = 1;
  opt.orderBatches.skullStripping = 2;
  
  % skull strip the mean image
  opt.skullStripMeanImg = 1; 
  
  % make matlab batch for segment and skullstip
  matlabbatch = [];
  matlabbatch = setBatchSegmentation(matlabbatch, opt, opt.funcMaskFileName);

  matlabbatch = setBatchSkullStripping(matlabbatch, BIDS, opt, subID);
  % run spm
  saveAndRunWorkflow(matlabbatch, 'meanImage_segment_skullstrip', opt, subID);


  % output the mask
  [meanImage, meanFuncDir] = getMeanFuncFilename(BIDS, subID, opt);
  %name the output accordingto the input image
  maskFileName = ['m' strrep(meanImage, '.nii', '_mask.nii')];
  mask = fullfile(meanFuncDir, maskFileName);



