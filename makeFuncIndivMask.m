% (C) Copyright 2020 RnB FFT-analysis developers

function mask = makeFuncIndivMask(opt)

  % function uses a mean functional image to create individual space mask

  % STEP 1
  % go to mricron and create skull stripped mean functional image
  % by using FSL BET function
  % in the future think about implementing FSL BET into matlab

  % read the mask image
  image = opt.funcMaskFileName;
  [imagePath, imageName, ext] = fileparts(image);

  % the output name for FSL BET image (skull stripped)
  betImageName = ['bet05_', imageName, ext];

  %SPM skull stripping - with Anat atm
  [~, opt, BIDS] = getData(opt);
  
  opt.orderBatches.segment = 1;
  opt.orderBatches.skullStripping = 2;
  
  % make matlab batch for segment and skullstip
  matlabbatch = [];
  matlabbatch = setBatchSegmentation(matlabbatch, opt, opt.funcMaskFileName);

  matlabbatch = setBatchSkullStripping(matlabbatch, BIDS, opt, subID);
  % run spm
  saveAndRunWorkflow(matlabbatch, 'meanImage_segment_skullstrip', opt, subID);


  % STEP 1.2
  % read already created bet05 image

  betImage = fullfile(imagePath, betImageName);

  % STEP 1.3
  % copy created bet05_meanfunc_mask images & rename with taskname
  if ~exist(betImage)

    rhythmBlockFuncDir = strrep(imagePath, opt.taskName, 'RhythmBlock');
    betOrigImage = strrep(fullfile(rhythmBlockFuncDir, betImageName), ...
                          opt.taskName, 'RhythmBlock');
    copyfile(betOrigImage, betImage);

  end

  % create mask name
  maskFileName = ['mask', betImageName];
  mask = fullfile(imagePath, maskFileName);

  % STEP 2
  if ~exist(mask)
    % Create a template & load the mask
    % A = load_untouch_nii('bet_05_meanuasub-pil001-PitchFT_run-001.nii');
    A = load_untouch_nii(betImage);

    C = A;
    C.fileprefix = 'C';
    C.img = [];

    idx = find(A.img > 0);
    A.img(idx) = 1;
    C.img = A.img;
    save_untouch_nii(C, mask);
  end

end

% function useBetFsl
%
%
% end
