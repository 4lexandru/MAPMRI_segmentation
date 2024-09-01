function kout = merge_seg(seg,kval)

% merge small voxels to neighboring clusters
% 
% Alexandru V. Avram, 
% NIH 2022
%

if nargin<2
    kval = 1e4; % dummy value 
end
% assign labels to spurious voxels denoted with 1e4
kout = seg;
disp(['Total of ',num2str(length(unique(seg(seg>0)))),' clusters found. Assigning small voxels to largest boundary neighbor...'])
[kfsz, kf] = bwconncomp_sizes(kout==kval);
for pp = 1:length(kfsz) % loop over small disconnected unlabelled clusters. assign new label based on largest boundary
    tmpz = 0*seg; tmpz(kf.PixelIdxList{pp}) = 1;
    aa = bwperim(1-tmpz); aa([1 end],:,:) = 0; aa(:,:,[1 end]) = 0; aa(:,[1 end],:) = 0;
    taa = kout(find(aa)); mv = max(taa(find(taa>0)));
    if ~isempty(mv)
        kout(kf.PixelIdxList{pp}) = mv;
    end
    disp([num2str(pp),' of ',num2str(length(kfsz))])
end
kout(kout==1e4) = 0;