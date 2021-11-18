function groupLevelzMapThreshold (opt, pvalue)
% find and read the group level zMap (output of FFT analysis), and
% thresholds it with a given pvalue, then saves it into the same folder
% with different name

%FWHM = opt.FWHM;

subNb = length(opt.subjects);

%cut off threshold for pvalue
threshold =  round(abs(norminv(pvalue)),2);
        

        % setup output directory
        fftDir = getFFTdir(opt, subLabel);
        [~, folder] = fileparts(fftDir);
        inputFolder = fullfile(fftDir, '..', 'group', folder);
        
        fftFileName = getZimage(inputFolder, subNb);

        % load image
        zHdr = spm_vol(fftFileName);
        zImg = spm_read_vols(zHdr);
        
        % threshold image with z-conversion of p<1e-4 
        zImg = zImg > threshold;
        
        % if option.save zmap correct save thresholded map
         % save zscore image
        if opt.save.zmap
            newName = [pvalue, fftFileName];
           saveThresholdedImage(zImg, zHdr, newName);
        end
        
        
end



function outputImage = getZimage(fftDir, subNb)
% find the zmap in a given folder

pattern = ['^whole-brain_AvgZTarget_subNb-', numb2str(subNb),'.nii$'];
outputImage = spm_select('FPList', fftDir, pattern);


end

function saveThresholdedImage (img, hdr, newName)

ext = '.nii';
newName = [newName,ext];

hdr.fname = spm_file(hdr.fname, 'filename', newName);

spm_write_vol(hdr, img);

end