function paths = getPaths()

    [ret, hostname] = system('hostname');

    paths = [];

    if strcmp(deblank(hostname), 'tux')
        paths.spm = '/home/tomo/Documents/MATLAB/spm12';
        paths.deriv = '/datadisk/data/RhythmCateg-fMRI/RhythmBlock';
    end
