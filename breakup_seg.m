function kout = breakup_seg(seg,th)

% splits up a segmentation into new unique disconnected labels
%
% Alexandru V. Avram, Ph.D.
% NIH 2022
%

if nargin<2
    th = 100; 
end

ct = 1; kout = seg*0;
CC = regionprops3(seg,'Centroid','VoxelIdxList','Volume');
 
for k = 1:length(CC.VoxelIdxList)
    tmp = seg*0; tmp(CC.VoxelIdxList{k}) = 1; bwconncomp(tmp,18);
    try 
        [consz, con] = bwconncomp_sizes(tmp); 
    disp(['---> ', num2str(length(find(consz>th))), ' regions above threshold out of ',num2str(length(consz))]);%
    for jj = 1:length(consz) % read all connected indices
        if consz(jj) > th
            kout(con.PixelIdxList{jj}) = ct; ct = ct+1;
        else
            kout(con.PixelIdxList{jj}) = 1e4; % assign a very large value placeholder for noise voxels
        end
    end
    end
end
disp(['Total of ',num2str(length(unique(kout(kout>0)))),' clusters found after initial spatial separation with threshold ',num2str(th),' voxels.'])

% combine all noise voxels and re-assign new labels
[kfsz, kf] = bwconncomp_sizes(kout==1e4);
disp(['Number of unlabelled voxels ', num2str(sum(kfsz)),'. Merging unlabelled small regions...'])
for pp = 1:length(kfsz) % loop over disconnected unlabelled clusters. assign new labels to clusters larger than th
    if kfsz(pp) > th
        kout(kf.PixelIdxList{pp}) = ct; ct = ct+1; %disp(ct);
    end
end
[kfsz, kf] = bwconncomp_sizes(kout==1e4);
disp(['Number of unlabelled voxels ', num2str(sum(kfsz)),'. Merging unlabelled small regions...'])