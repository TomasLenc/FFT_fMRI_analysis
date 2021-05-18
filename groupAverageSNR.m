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
    
    
    % get Target nii files
    % whole-brain_AvgZTarget_s2wuasub-001_ses-001_task-PitchFT_bold.nii
    avgZFileName = [maskType, '_AvgZTarget_', boldFileName, '.nii'];
    avgZFileFolder{iSub} = fullfile(inputDir,avgZFileName);
    
    % get Ratio nii
    % whole-brainAvgRatioTarget_s2wuasub-001_ses-001_task-PitchFT_bold
    
    
    % get .mat FT files
    % whole-brain_AvgZTarget_s2wuasub-001_ses-001_task-PitchFT_bold_FT.mat
    ftFileName = [avgZFileName(1:end-4),'_FT', '.mat'];
    ftFileFolder{iSub} = fullfile(inputDir, ftFileName);
    
end

disp('Nii Files:');
for iFile = 1:length(avgZFileFolder)
    disp([avgZFileFolder{iFile}]);
end

fprintf(' \n NumSubjects: %i  \n\n', numel(opt.subjects));

disp('FT Files:');
for iFile = 1:length(ftFileFolder)
    disp([ftFileFolder{iFile}]);
end

% Empty matrix of 4 dimensions (first 3 dimensions are the brain image,
% the fourth dimention is the subject number)
z = [];

% first subject Number
iSub = 1;

% loop for each subject
while iSub <= numel(opt.subjects)
    
    % load the searchlight map
    hdr = spm_vol(avgZFileFolder{iSub});
    img = spm_read_vols(hdr);
    
    fprintf('Loading of Map %.0f finished. \n', iSub);
    
    % concatenate each subject to the 4th dimension
    z = cat(4, z, img);
    
    FT = load(ftFileFolder{iSub});
    
    % increase the counter
    iSub = iSub + 1;
end

%% Mean Accuracy
% calcuate mean accuracy maps
% Calculate mean of each voxel across subjects (4th dimension)
means = [];
means = mean(z, 4);
adjustMeanImg = means ./(sqrt(numel(opt.subjects)));
meanImg = means;

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
adjustMeanHdr = hdr;
newFileName = [maskType, '_AdjustAvgZTarget_subNb-', num2str(numel(opt.subjects)), '.nii'];
adjustMeanHdr.fname = spm_file(adjustMeanHdr.fname, 'filename', newFileName);

% save result as .nii file
spm_write_vol(adjustMeanHdr, adjustMeanImg);


% %% setup parameters for FFT analysis
% % mri.repetition time(TR)
% repetitionTime = 1.75;
% 
% % repetition of steps/categA
% patternDuration     = 12 * 0.190;
% segmentDuration     = 4 * patternDuration;
% stepDuration        = opt.nStepsPerPeriod * segmentDuration;
% 
% % Number of vol before/after the rhythmic sequence (exp) are presented
% onsetDelay = 2;
% endDelay = 4;
% 
% % use neighbouring 4 bins as noise frequencies
% cfg.binSize = 4;
% cfg.gap = 1;
% 
% N = 104; % after resampling
% 
% % calculate frequencies
% oddballFreq = 1 / stepDuration;
% oldFs =  1 / repetitionTime;
% fs = 1 / (182.4 / N);
% 
% % frequencies
% freq = fs / 2 * linspace(0, 1, N / 2 + 1);
% 
% % target frequency (this is a bad name, because it's a freq. bin INDEX)
% cfg.targetFrequency = round(N * oddballFreq / fs + 1);
% 
% % harmonics of the target frequency and their bin indices
% cfg.harmonics = freq(cfg.targetFrequency) * opt.whichHarmonics;
% cfg.idxHarmonics = dsearchn(freq', cfg.harmonics');
% 
% % number of bins for phase histogram
% cfg.histBin = 20;
% 
% % threshold for choosing voxels for the phase distribution analysis
% cfg.thresh = 4;

end
