function mXavgHarmonics = getAvgHarmonics(mX, idxHarmonics, binSize)

  % if 1D make it a column vector
  if size(mX, 1) == 1
    mX = mX';
  end

  nHarmonics = length(idxHarmonics);
  nVox = size(mX, 2);

  mXavgHarmonics = nan(binSize * 2 + 1, length(nHarmonics), nVox);

  for iHarm = 1:nHarmonics
    idxMin = idxHarmonics(iHarm) - binSize;
    idxMax = idxHarmonics(iHarm) + binSize;

    mXavgHarmonics(:, iHarm, :) = mX([idxMin:idxMax], :);

  end

  mXavgHarmonics = squeeze(mean(mXavgHarmonics, 2));
