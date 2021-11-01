function [opt] = getMaskFile(opt)

% this is a small function to get mask according to the input.
% either whole-brain masks will be uploaded or  freesurfer, hmat binary
% maps. 
% neurosythn maps should need extra care - not sure if it is working

% according to the maskType, certain masks will be loaded
maskType = opt.maskType;

%preallocate
maskFile = cell(numel(opt.subjects), 1);

if strcmp(maskType, 'whole-brain')
    % create a whole brain functional mean image mask
    % so the mask will be in the same resolution/space as the functional images
    % one may not need it if they are running bidsFFX since it creates a
    % mask.nii by default
    % opt.skullstrip.threshold = 0.5; -->provides bigger mask thn default value
    opt.skullstrip.mean = 1;
    
    maskFile = bidsWholeBrainFuncMask(opt);
    
    opt.maskLabel = maskType;
    
else
    
    % find the rois
    opt = chooseMask(opt, maskType);
    
    % assign ROI folder
    maskPath = opt.maskPath;
    maskName = opt.maskName;   
    
    if strcmp(opt.space, 'individual')
        
        for iSub = 1:numel(opt.subjects)
            
            % read subject ID
            subFolder = ['sub-', opt.subjects{iSub}];
            
            for iMask = 1:numel(maskName)
                maskFile{iSub, iMask}  = fullfile(opt.maskPath, ...
                    subFolder, ...
                    [opt.maskBaseName, ...
                    opt.maskName{iMask}]);
            end
            
        end
        
    elseif strcmp(opt.space, 'MNI')
        
        for iSub = 1:numel(opt.subjects)
            for iMask = 1:numel(maskName)
                maskFile{iSub, iMask} = fullfile(maskPath, maskName{iMask});
            end
        end
        
    end
    
    
    
end

% assign/save the mask files into structure
opt.maskFile = maskFile;

end
