function [coord] = getVoxelCoordinate(hdr, img, maskImg, voxelNbToPlot)
    % gets the best N voxels coordinate by sorting the img values in descending
    % order

    % img is the 3D matrix (image)
    % hdr : header of the img
    % number of top N voxels : voxelNbToPlot

    % an example call for the vx coordinate:
    % coord.voxelSpaceXyz(iVox,:) = [25 45 39]

    % Gets the Voxel space to World space transformation matrix of this image
    transformationMatrix = hdr(1).mat;

    % sort the z-value 1D map to find max N voxels to plot
    [valuesSorted, idxSorted] = sort(img(:), 'descend');

    %  control of NaN/zeros
    idxNoNans = ~isnan(valuesSorted);
    valuesSorted = valuesSorted(idxNoNans);
    idxSorted = idxSorted(idxNoNans);

    for iVox = 1:voxelNbToPlot

        % convert linear indices into 3D indices
        [x, y, z] = ind2sub(size(img), idxSorted(iVox));

        % have to take the transpose (') of the voxel indices vector
        transposeVoxelSubscripts = [x; y; z];

        % pad it with an extra one to be able to multiply it with the tranformation matrix.
        worldSpaceXyz = transformationMatrix * [transposeVoxelSubscripts; 1];

        % Only the three first value of this vector are of interest to us
        worldSpaceXyz = worldSpaceXyz(1:3);

        %         %% an example from vx script to MNI space
        %         transformationMatrix = boldHdr(1).mat;
        %         MNIcoord = cor2mni([x,y,z], transformationMatrix);

        % save into a struct
        coord.voxelSpaceXyz(iVox, 1) = x;
        coord.voxelSpaceXyz(iVox, 2) = y;
        coord.voxelSpaceXyz(iVox, 3) = z;
        coord.index(iVox) = idxSorted(iVox);
        coord.worldSpaceXyz(iVox, 1) = worldSpaceXyz(1);
        coord.worldSpaceXyz(iVox, 2) = worldSpaceXyz(2);
        coord.worldSpaceXyz(iVox, 3) = worldSpaceXyz(3);
        coord.value(iVox) = valuesSorted(iVox);

    end

    % convert to 1D linearized indices on the masked image
    coord.indexMasked = img3DtoMask1D(coord.voxelSpaceXyz, maskImg);

    % an example is below:

    % %% Getting the world space coordinate of a given voxel
    % VoxelSubscripts = [x y z];
    %
    % % Gets the Voxel space to World space transformation matrix of this image
    % transformationMatrix = boldHdr(1).mat;
    %
    % % The dimension of each voxel are the 3 first values of the main diagonal of this matrix
    % % Get the value along the main diagonal of the matrix and stores it in a
    % % temporary variable
    % temp = diag(transformationMatrix);
    %
    % % We have to take the transpose (') of the voxel indices vector and pad it with an
    % % extra one to be able to multiply it with the tranformation matrix.
    % temp = [VoxelSubscripts' ; 1];
    % worldSpaceXyz = transformationMatrix * temp;
    %
    % % Only the three first value of this vector are of interest to us
    % worldSpaceXyz = worldSpaceXyz(1:3);

end

function [idxImg3D, idxImg1D] = mask1DtoImg3D(idxMask1D, maskImg)
    % convert from linearized 1D index of a voxel in masked data to
    % 3D coordinate index in the original image (before mask was applied)
    idxImg3D = nan(length(idxMask1D), 3);
    idxImg1D = nan(1, length(idxMask1D));

    for i = 1:length(idxMask1D)
        % make a dummy vector of the linearized 1D mask data
        maskData = zeros(1, sum(maskImg(:) > 0));

        % put the requested index to the corresponding positon
        maskData(idxMask1D(i)) = 1;

        % get dimensions & allocate 3-D img
        img = zeros(size(maskImg));

        img(find(maskImg > 0)) = maskData;

        % we can also get the 1D index in the original image if we need
        idxImg1D(i) = find(img);

        % find the 3D coordinates x,y,z
        [x, y, z] = ind2sub(size(maskImg), idxImg1D(i));
        idxImg3D(i, :) = [x, y, z];
    end
end

function idxMask1D = img3DtoMask1D(idxImg3D, maskImg)
    % convert from 3D coordinate index of the original image to linearized 1D
    % index in a masked image
    idxMask1D = nan(1, size(idxImg3D, 1));

    for i = 1:size(idxImg3D, 1)
        img = zeros(size(maskImg));

        idxImg1D = sub2ind(size(maskImg), idxImg3D(i, 1), idxImg3D(i, 2), idxImg3D(i, 3));
        img(idxImg1D) = 1;

        maskedImg = img(maskImg > 0);

        idxMask1D(i) = find(maskedImg);
    end
end
