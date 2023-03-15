function f = plotmXBestVox(freq, mX, coord, idxHarmonics,  varargin)
    % This function finds voxels with highest targetVal and plots ther magnitude
    % spectra (mX) as a function of frequency (freq)
    %
    % Args:
    % -----
    %     freq:       [1xN] array with frequency value in Hz for each frequency bin
    %                 in the spectra
    %     mX:         [frequency x voxel] array of spectral magnitudes
    %     coord :     ...
    %     idxHarmonics : index, which frequency bins correspont to target frequency
    %                 (and harmonics if relevant)
    %     varargin:
    %                 'ratio' : flag if y-scale should be adjusted around 1 for ratio values
    %
    %
    % Returns:
    % --------
    %     f : figure handle

    % if 1D make it a column vector
    if size(mX, 1) == 1
        mX = mX';
    end

    % get number of best voxeld we'll be plotting, each one will be in a
    % separate subplot
    nVox = length(coord.index);

    % make figure, pack subplots
    f = figure('color', 'white', 'Position', [123 782 555 800]);
    pnl = panel(f);
    pnl.pack('v', nVox);

    for iVox = 1:nVox

        % Get the magnitude spectra for the current voxel. Use the linear 1D
        % index for the masked data. We need this because we don't have FFT for
        % all the voxels in the image, we only calculated it on the masked
        % voxels...
        mXbest = mX(:, coord.indexMasked(iVox));
        mXbest(1) = 0;
        mXbest = mXbest(1:length(freq));

        pnl(iVox).select();

        % plot all the magnitude spectra in grey
        stem(freq, mXbest, ...
             'marker', 'none', ...
             'color', [0.6, 0.6, 0.6], ...
             'linew', 2);

        hold on;
        % highlight the target frequency in red
        stem(freq(idxHarmonics(1)), mXbest(idxHarmonics(1)), ...
             'marker', 'none', ...
             'color', 'r', ...
             'linew', 2);

        % highlight targett frequency harmonics in organge
        stem(freq(idxHarmonics(2:end)), mXbest(idxHarmonics(2:end)), ...
             'marker', 'none', ...
             'color', [255, 138, 102] / 255, ...
             'linew', 2);

        % get the y-value at target frequency (and harmonics)
        vals = mXbest(idxHarmonics);
        % select the maximum y-value to set limits
        maxVal = max(vals);
        % based on the max y-value of interest, set reasonable y limits
        % (we do this separately for each subplot)
        if any(strcmpi(varargin, 'ratio'))
            ylim([1, 1.1 * maxVal]);
            valueLabel = 'snr';
        else
            ylim([-0.1 * maxVal, 1.1 * maxVal]);
            valueLabel = 'z';
        end

        % make it pretty
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.XTick = [0.2];
        ax.YAxis.TickDirection = 'out';
        ax.XAxis.TickDirection = 'out';
        if iVox < nVox
            ax.XTickLabel = {};
        end

        title(sprintf('%s=%.2f  vox=[%d %d %d]', ...
                      valueLabel, ...
                      coord.value(iVox), ...
                      coord.voxelSpaceXyz(iVox, :)));
    end

    pnl.de.marginbottom = 10;
    pnl.de.margintop = 1;

    pnl.marginbottom = 7;
