function mXsnr = baselineCorrect(mX, cfg, varargin)

    if any(strcmpi(varargin, 'fun'))
        blfun = varargin{ find(strcmpi(varargin, 'fun')) + 1 };
    else
        blfun = @minus;
    end

    gap = cfg.gap;
    binSize = cfg.binSize;

    % if 1D make it a column vector
    if size(mX, 1) == 1
        mX = mX';
    end

    mXsnr = zeros(size(mX));
    N = size(mX, 1);

    for i = 1:N
        idx11 = max(1, i - gap - binSize / 2);
        idx12 = max(1, i - gap - 1);
        idx21 = min(N, i + gap + 1);
        idx22 = min(N, i + gap + binSize / 2);

        signal = mX(i, :);

        noiseMean = mean(mX([idx11:idx12, idx21:idx22], :), 1);

        mXsnr(i, :) = blfun(signal, noiseMean);
    end
