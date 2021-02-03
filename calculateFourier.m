function [targetSNR, cfg] = calculateFourier(X, Xraw, cfg)
  % Fourier analysis of fMRI time series data, returns the SNR at a given
  %        frequency for each voxel
  %
  % Xiaoqing Gao, Dec 6, 2017, Louvain-la-Neuve, dr.x.gao@gmail.com
  %
  % Input data:
  %        1.X: rows are MRI measurements (at time points); columns are
  %          voxels; Xraw, raw data, no linear detrending

  %        2.TargetFrequency: the frequency bin (ordinal number) of the
  %          target frequency

  %        3.BinSize: number of frequency bins surronding the target
  %          frequency bin (e.g., 40 means ?20 bins, with a gap of 1 bin from
  %          the target frequncy). These neighbouring frequencies are treaed as
  %          noise.

  %        4.Thresh: a threshold value (e.g., 3.719). The distribution of
  %          the phase values of the voxels above this threshold value are
  %          used to define activation/deactivation.

  %        5.histBin: initial number of bins used for calculating the histogram
  %          of phase distribution for defining activation/deactivation, if
  %          there are more than one bins with maximum count, histBin+1
  %
  % Output: 1. TargetSNR: SNR (z-score) at the target frequency for each
  %           voxels

  %        2. TargetPhase: Fourier phase of the target frequency for each
  %           voxel

  %        3. TargetSNRsigned: Activation (+)/deactivation (-) defined by
  %           phase values and applied to the SNR values for each voxel. This
  %           can be used to generate a signed map.

  %        4. tSNR

  % Steps of the analysis
  % 1. FFT of the time series

  targetFrequency = cfg.targetFrequency;
  binSize = cfg.binSize;
  thresh = cfg.thresh;
  histBin = cfg.histBin;

  tSNR = mean(Xraw) ./ std(Xraw);

  FT = fft(X);

  % 2. define noise frequencies based on TargetFrequency and BinSize with a
  % gap of 1
  gap = 1;
  noiseFs = [(targetFrequency - binSize / 2 - gap): ...
             (targetFrequency - 1 - gap) (targetFrequency + 1 + gap): ...
             (targetFrequency + binSize / 2 + gap)];

  % 3. calculate the mean and SD of the amplitudes of the noise frequencies
  FTNoise = FT(noiseFs, :);
  AmpNoise = abs(FTNoise);
  NoiseMean = mean(AmpNoise, 1);
  NoiseSD = std(AmpNoise, 0, 1);

  % 4. calculate SNR (z-score) of the target frequency based on the mean and SD of the
  % noise frequencies
  targetSNR = (abs(FT(targetFrequency, :)) - NoiseMean) ./ NoiseSD;

  % 5. using the distribution of phase of the target frequency to define the sign
  targetPhase = angle(FT(targetFrequency, :));

  % 5.1 find the peak of the distribution of positive phase
  % It assums that in the experiment, stimulus onset is at 0 phase
  TargetPhaseP = targetPhase(targetPhase > 0); % positive phase values

  while 1
    [n, x] = hist(TargetPhaseP(targetSNR(targetPhase > 0) > thresh), ...
                  histBin);
    xcenter = x(n == max(n));
    if length(xcenter) > 1
      histBin = histBin + 1;
    else
      break
    end
  end

  % 5.2 assign the peak phase ? pi/2 to be 1 and the others to be -1
  phaseDiff = abs(targetPhase - xcenter);
  phaseIndex = zeros(size(phaseDiff));
  phaseIndex(phaseDiff <= (pi / 2)) = 1;
  phaseIndex(phaseDiff > (pi / 2)) = -1;

  % 5.3 apply sign to target SNR
  targetSNRsigned = targetSNR .* phaseIndex;

  % unused parameters for now
  cfg.targetPhase = targetPhase;
  cfg.targetSNRsigned = targetSNRsigned;
  cfg.tSNR = tSNR;
  %
