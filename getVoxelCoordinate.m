function [coord] = getVoxelCoordinate(hdr, img, voxelNbToPlot)
% gets the best N voxels coordinate by sorting the img values in descending
% order

% img is the 3D matrix (image)
% hdr : header of the img
% number of top N voxels : voxelNbToPlot

% Gets the Voxel space to World space transformation matrix of this image
transformationMatrix = hdr(1).mat;

% sort the z-value 1D map to find max N voxels to plot
[zValuesSorted, idxSorted] = sort(img(:), 'descend');

for iVox = 1:voxelNbToPlot
    
    % convert linear indices into 3D indices
    [x, y, z] = ind2sub(size(img),idxSorted(iVox));
    
    % have to take the transpose (') of the voxel indices vector
    transposeVoxelSubscripts = [x; y; z];
    
    % pad it with an extra one to be able to multiply it with the tranformation matrix.
    worldSpaceXyz = transformationMatrix * [transposeVoxelSubscripts ; 1];
    
    % Only the three first value of this vector are of interest to us
    worldSpaceXyz = worldSpaceXyz(1:3);
    
    %         %% an example from vx script to MNI space
    %         transformationMatrix = boldHdr(1).mat;
    %         MNIcoord = cor2mni([x,y,z], transformationMatrix);
    
    % save into a struct
    coord.voxelSpaceXyz(iVox,1) = x;
    coord.voxelSpaceXyz(iVox,2) = y;
    coord.voxelSpaceXyz(iVox,3) = z;
    coord.index(iVox) = idxSorted(iVox);
    coord.worldSpaceXyz(iVox,1) = worldSpaceXyz(1);
    coord.worldSpaceXyz(iVox,2) = worldSpaceXyz(2);
    coord.worldSpaceXyz(iVox,3) = worldSpaceXyz(3);
    coord.zValue(iVox) = zValuesSorted(iVox);
    
    % to call the vx coordinate:
    % coord.voxelSpaceXyz(iVox,:) = [25 45 39]
    
end

%         % vector array structure
%             coord(iVox).voxelSpace.x = x;
%             coord(iVox).voxelSpace.y = y;
%             coord(iVox).voxelSpace.z = z;
%             coord(iVox).index = idxSorted(iVox);
%             coord(iVox).worldSpace.x = worldSpaceXyz(1);
%             coord(iVox).worldSpace.y = worldSpaceXyz(2);
%             coord(iVox).worldSpace.z = worldSpaceXyz(3);
%             coord(iVox).zValue = zValues(iVox);


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
% 
% coord(iVox).voxelSpace.x = x;
% coord(iVox).voxelSpace.y = y;
% coord(iVox).voxelSpace.z = z;
% coord(iVox).index = idxSorted(iVox);
% coord(iVox).worldSpace.x = worldSpaceXyz(1);
% coord(iVox).worldSpace.y = worldSpaceXyz(2);
% coord(iVox).worldSpace.z = worldSpaceXyz(3);
% coord(iVox).zValue = zValuesSorted(iVox);
