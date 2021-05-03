function f = plotmXBestVox(freq, mX, targetVal, nPlot, idxHarmonics, coord, varargin)

    coord.voxelSpaceXyz 
    
    % if 1D make it a column vector
    if size(mX, 1) == 1
        mX = mX';
    end

    [~, idxSorted] = sort(targetVal, 'descend');

    idxSorted(isnan(mX(idxSorted))) = [];

    f = figure('color', 'white', 'Position', [123 782 555 800]);
    pnl = panel(f);
    pnl.pack('v', nPlot);

    for iVox = 1:nPlot

        mXbest = mX(:, idxSorted(iVox));
        mXbest(1) = 0;
        mXbest = mXbest(1:length(freq));

        pnl(iVox).select();

        stem(freq, mXbest, ...
             'marker', 'none', ...
             'color', [0.6, 0.6, 0.6], ...
             'linew', 2);

        hold on;
        stem(freq(idxHarmonics(1)), mXbest(idxHarmonics(1)), ...
             'marker', 'none', ...
             'color', 'r', ...
             'linew', 2);

        stem(freq(idxHarmonics(2:end)), mXbest(idxHarmonics(2:end)), ...
             'marker', 'none', ...
             'color', [255, 138, 102] / 255, ...
             'linew', 2);

        vals = mXbest(idxHarmonics);
        maxVal = max(vals);

        if any(strcmpi(varargin, 'ratio'))
            ylim([1, 1.1 * maxVal]);
        else
            ylim([-0.1 * maxVal, 1.1 * maxVal]);
        end

        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.XTick = [0.2];
        ax.YAxis.TickDirection = 'out';
        ax.XAxis.TickDirection = 'out';
        if iVox < nPlot
            ax.XTickLabel = {};
        end

        title(sprintf('voxel %d', idxSorted(iVox)));
    end

    pnl.de.marginbottom = 10;
    pnl.de.margintop = 1;

    pnl.marginbottom = 7;
