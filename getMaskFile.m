function [maskFile, maskLabel] = getMaskFile(opt)

    maskType = opt.maskType;

    maskFile = cell(numel(opt.subjects), 1);

    if strcmp(maskType, 'whole-brain')
        % create a whole brain functional mean image mask
        % so the mask will be in the same resolution/space as the functional images
        % one may not need it if they are running bidsFFX since it creates a
        % mask.nii by default
        % opt.skullstrip.threshold = 0.5; -->provides bigger mask thn default value
        opt.skullstrip.mean = 1;
        maskFile = bidsWholeBrainFuncMask(opt);

        maskLabel = maskType;
    end

    for iSub = 1:numel(opt.subjects)

        % read subject ID
        subLabel = opt.subjects{iSub};

        % freesurfer or neurosynth ROIs?
        if strcmp(maskType, 'freesurfer') && strcmp(opt.space, 'individual')

            % check which ROI folder
            maskPath = fulllfile(opt.roiDir, ...
                                 [maskType, '/sub-', subLabel]);

            maskName = {'thres5_s1_dec_rlauditorycx.nii', ...
                        'thres5_s1_dec_rrauditorycx.nii', ...
                        'rlbasalganglia.nii', ...
                        'rrbasalganglia.nii'};

            maskLabel = {'leftAud', 'rightAud', 'leftBG', 'rightBG'};

            for iMask = 1:numel(maskName)

                maskFile{iSub, iMask} = fullfile(maskPath, maskName{iMask});
            end

        elseif strcmp(maskType, 'neurosynth') && strcmp(opt.space, 'MNI')

            maskPath = fullfile(opt.roiDir, maskType, 'functional', ...
                                'derivatives');

            maskName = {'leftbin_rnativeThres_7_auditory_FDR_0.01.nii', ...
                        'rightbin_rnativeThres_6_auditory_FDR_0.01.nii', ...
                        'rrthres_7sma_FDR_0.01.nii', ...
                        'leftrrthres_5premotor_FDR_0.01.nii', ...
                        'rightbin_rnativeThres_5_premotor_FDR_0.01.nii'};

            % use in output labeling
            maskLabel = {'leftAud', 'rightAud', 'SMA', 'leftPre', 'rightPre'};

            for iMask = 1:numel(maskName)
                maskFile{iSub, iMask} = fullfile(maskPath, maskName{iMask});
            end

        end

    end

    %% find the "biggest" ROIs to improvie the coverage of the mask

    %   % use in output name
    %   roiSource = 'neurosnyth';
    %
    %   maskPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
    %                       '..', '..', '..', 'RhythmCateg_ROI', 'neurosynth', ...
    %                       'functional', 'derivatives');
    %
    %   % masks to decode
    %   % maskName = {'leftrrthres_7premotor_FDR_0.01.nii', ...
    %   %            'rightrrthres_7premotor_FDR_0.01.nii',...
    %   %            'rrthres_10sma_FDR_0.01.nii'};
    %   maskName = {'leftbin_rnativeThres_7_auditory_FDR_0.01.nii', ...
    %               'rightbin_rnativeThres_6_auditory_FDR_0.01.nii', ...
    %               'rrthres_7sma_FDR_0.01.nii', ...
    %               'leftrrthres_5premotor_FDR_0.01.nii', ...
    %               'rightbin_rnativeThres_5_premotor_FDR_0.01.nii'};
    %
    %   % use in output roi name
    %   maskLabel = {'leftAud', 'rightAud', 'SMA', 'leftPremotor', 'rightPremotor'};
    %   % maskLabel = {'leftPremotor','rightPremotor', 'SMA'};
    %
    %   % parcels
    %   % if use parcels, re-writes mask names:
    %   if opt.mvpa.useParcel == 1
    %
    %     roiSource = 'freesurfer';
    %
    %     maskPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
    %                         '..', '..', '..', 'RhythmCateg_ROI', 'freesurfer');
    %
    %     maskName = {'thres5_s1_dec_rlauditorycx.nii', ...
    %                 'thres5_s1_dec_rrauditorycx.nii', ...
    %                 'rlbasalganglia.nii', ...
    %                 'rrbasalganglia.nii'};
    %
    %     maskLabel = {'leftAud', 'rightAud', 'leftBG', 'rightBG'};

end
