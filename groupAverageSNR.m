function opt = groupAverageSNR(opt)

%% set up experiment related info

% number of steps per analysed period
% Set default for RhythmFT and PitchFT analyses
opt.whichHarmonics = [1, 2];

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
end

% setup output directory
fftDir = fullfile(opt.derivativesDir, '..', 'rnb_fft');
destinationDir = fullfile(fftDir, 'group');

% get mask image
% use a predefined mask, only calculate voxels within the mask
% below is same resolution as the functional images
maskType = opt.maskType;

%% let's start
for iSub = 1:numel(opt.subjects)
    
    %get subject label
    subLabel = opt.subjects{iSub};
    
    % input directory
    inputDir = createOutputDirectory(opt, subLabel);
    [~,folder] = fileparts(inputDir);
    
    % input midfile name
    opt = getSpecificBoldFiles(opt, subLabel);
    [~, boldFileName, ~] = fileparts(opt.allFiles{1});
    boldFileName = regexprep(boldFileName, 'run-(\d*)_', '');
    
    
    
if strcmpi(opt.taskName, 'RhythmBlock')
    
    % get Target nii files
    avgZFileName = ['AvgZTarget_', boldFileName, '.nii'];
    
    % get ratio target
    ratioFileName = ['AvgRatioTarget_', boldFileName, '.nii'];
else
    
    % get Target nii files
    avgZFileName = [maskType, '_AvgZTarget_', boldFileName, '.nii'];

    % get ratio target
    ratioFileName = [maskType, 'AvgRatioTarget_', boldFileName, '.nii'];
    
end

% keep these names/.nii files
avgZFileFolder{iSub} = fullfile(inputDir,avgZFileName);
ratioFileFolder{iSub} = fullfile(inputDir,ratioFileName);

end

disp('Nii Files:');
for iFile = 1:length(avgZFileFolder)
    disp([avgZFileFolder{iFile}]);
end

fprintf(' \n NumSubjects: %i  \n\n', numel(opt.subjects));

disp('Nii Files:');
for iFile = 1:length(ratioFileFolder)
    disp([ratioFileFolder{iFile}]);
end


% Empty matrix of 4 dimensions (first 3 dimensions are the brain image,
% the fourth dimention is the subject number)
z = [];
zratio = [];

% first subject Number
iSub = 1;

% loop for each subject
while iSub <= numel(opt.subjects)
    
    % load the average z-score of target frequency
    hdr = spm_vol(avgZFileFolder{iSub});
    img = spm_read_vols(hdr);
    
    fprintf('Loading of Map %.0f finished. \n', iSub);

    % concatenate each subject to the 4th dimension
    z = cat(4, z, img);
    
    % load the average ratio
    hdrRatio = spm_vol(ratioFileFolder{iSub});
    imgRatio = spm_read_vols(hdrRatio);
    
    fprintf('Loading of Map %.0f finished. \n', iSub);

    % concatenate each subject to the 4th dimension
    y = cat(4, y, imgRatio);

    % increase the counter
    iSub = iSub + 1;
end

%% Mean Accuracy
% calcuate mean accuracy maps
% Calculate mean of each voxel across subjects (4th dimension)
means = [];
means = mean(z, 4);
meansRatio = mean(y, 4);

adjustMeanImg = means ./(sqrt(numel(opt.subjects)));
meanImg = means;
meanRatioImg = meansRatio;

meanHdr = hdr;


newFolderToSave = fullfile(destinationDir,folder);
if ~exist(newFolderToSave,'dir')
    mkdir(newFolderToSave);
end

newFileName = [maskType, '_AvgZTarget_subNb-', num2str(numel(opt.subjects)), '.nii'];

meanHdr.fname = spm_file(meanHdr.fname, 'path', newFolderToSave);
meanHdr.fname = spm_file(meanHdr.fname, 'filename', newFileName);

% save result as .nii file
spm_write_vol(meanHdr, meanImg);


% adjusted mean map
% The averaged z-scores are no longer from a standard normal distribution. 
% Instead, they have a standard deviation of 1??n (n = numSubjects).
adjustMeanHdr = meanHdr;
newFileName = [maskType, '_AdjustAvgZTarget_subNb-', num2str(numel(opt.subjects)), '.nii'];
adjustMeanHdr.fname = spm_file(adjustMeanHdr.fname, 'filename', newFileName);

% save result as .nii file
spm_write_vol(adjustMeanHdr, adjustMeanImg);


% save ratio
meanRatioHdr = meanHdr;
newFileName = [maskType, '_AvgRatio_subNb-', num2str(numel(opt.subjects)), '.nii'];
adjustMeanHdr.fname = spm_file(adjustMeanHdr.fname, 'filename', newFileName);

% save result as .nii file
spm_write_vol(meanRatioHdr, meanRatioImg);

end
