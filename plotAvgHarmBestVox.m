function f = plotAvgHarmBestVox(mXavgHarmonics, coord)

    % get number of best voxeld we'll be plotting, each one will be in a
    % separate subplot
    nVox = length(coord.index);

    f = figure('color', 'white', 'Position', [131 728 1744 140]);
    pnl = panel(f);
    pnl.pack('h', nVox);

    for iVox = 1:nVox

        mXbest = mXavgHarmonics(:, coord.indexMasked(iVox));

        pnl(iVox).select();

        h = stem(mXbest, ...
                 'marker', 'none', ...
                 'color', [0.6, 0.6, 0.6], ...
                 'linew', 4);
        hold on;

        centerIdx = floor(length(mXbest) / 2) + 1;

        h = stem(centerIdx, mXbest(centerIdx), ...
                 'marker', 'none', ...
                 'color', 'r', ...
                 'linew', 4);

        title(sprintf('z=%.2f  vox=[%d %d %d]', ...
                      coord.value(iVox), ...
                      coord.voxelSpaceXyz(iVox, :)));

        set(gca, 'xtick', []);
    end

    pnl.ylabel('magnitude');
    pnl.marginbottom = 3;
