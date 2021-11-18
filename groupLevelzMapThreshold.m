function groupLevelzMapThreshold (opt, pvalue)
% find and read the group level zMap (output of FFT analysis), and
% thresholds it with a given pvalue, then saves it into the same folder
% with different name

%FWHM = opt.FWHM;

subNb = length(opt.subjects);

%cut off threshold for pvalue
threshold =  round(abs(norminv(pvalue)),2);

% prepare pvalue for labeling the output
pvalue = strrep(num2str(pvalue), '.', '');

% dummy call
subLabel = opt.subjects{1};

% setup output directory
fftDir = getFFTdir(opt, subLabel);
[~, folder] = fileparts(fftDir);
inputFolder = fullfile(fftDir, '..', '..', '..','group', folder);

fftFile = getZimage(inputFolder, subNb);

% load image
zHdr = spm_vol(fftFile);
zImg = spm_read_vols(zHdr);

% threshold image with z-conversion of p<1e-4
% zImg = zImg > threshold; % this is binary conversion
idx = zImg < threshold;
zImg(idx) = 0;

% if option.save zmap correct save thresholded map
% save zscore image
if opt.save.zmap
    [~, name , ~] = fileparts(fftFile);
    newName = [name, '_thres-p', pvalue];
    saveThresholdedImage(zImg, zHdr, newName);
end


end



function outputImage = getZimage(fftDir, subNb)
% find the zmap in a given folder

pattern = ['^whole-brain_AvgZTarget_subNb-', num2str(subNb),'.nii$'];
outputImage = spm_select('FPList', fftDir, pattern);


end

function saveThresholdedImage (img, hdr, newName)

ext = '.nii';
newName = [newName,ext];

hdr.fname = spm_file(hdr.fname, 'filename', newName);

spm_write_vol(hdr, img);

end