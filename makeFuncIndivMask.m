% (C) Copyright 2020 RnB FFT-analysis developers

function mask = makeFuncIndivMask(opt)

    % function cpp-spm repo/functions and create skull stripped
    % mean functional image
    
    % read the dataset
    [~, opt, BIDS] = getData(opt);
    subID = opt.subjects{1};
    
    % call/create the mask name
    [meanImage, meanFuncDir] = getMeanFuncFilename(BIDS, subID, opt);
    
    % name the output accordingto the input image
    maskFileName = ['m' strrep(meanImage, '.nii', '_mask.nii')];
    mask = fullfile(meanFuncDir, maskFileName);

    % ask if mask exist, if not create it:
    if ~exist(mask,'file')
        
        % set batch order since there is dependencies
        opt.orderBatches.segment = 1;
        opt.orderBatches.skullStripping = 2;

        % opt for running skull strip on the mean image
        opt.skullStripMeanImg = 1;

        % make matlab batch for segment and skullstip
        matlabbatch = [];
        matlabbatch = setBatchSegmentation(matlabbatch, opt, opt.funcMaskFileName);

        matlabbatch = setBatchSkullStripping(matlabbatch, BIDS, opt, subID);
        % run spm
        saveAndRunWorkflow(matlabbatch, 'meanImage_segment_skullstrip', opt, subID);
    end
