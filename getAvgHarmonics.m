function mXavgHarmonics = getAvgHarmonics(mX, idxHarmonics, binSize)
% This function takes the whole spectra and cuts out segments around each
% frequency (harmonic), and averages these segments. 
% 
% Arguments: 
% ----------
%     mX :            array [frequency x voxel] of magnitude spectra
%     idxHarmonics :  indices specifying bins of frequencies of interest 
%     binSize :       number of bins on either size of each frequency of interest
%                     that will be segmented
% 
% Returns: 
% --------
%     mXavgHarmonics : array [voxel x bin], segment averaged overt frequencies 
%                      of interest
% ========================================================================

% if 1D make it a column vector
if size(mX, 1) == 1
    mX = mX';
end

nHarmonics = length(idxHarmonics);
nVox = size(mX, 2);

% allocate 
mXavgHarmonics = nan(binSize * 2 + 1, length(nHarmonics), nVox);

% go over frequencies of interest 
for iHarm = 1:nHarmonics
    idxMin = idxHarmonics(iHarm) - binSize;
    idxMax = idxHarmonics(iHarm) + binSize;

    mXavgHarmonics(:, iHarm, :) = mX([idxMin:idxMax], :);

end

% average over freqencies of interest 
mXavgHarmonics = squeeze(mean(mXavgHarmonics, 2));
