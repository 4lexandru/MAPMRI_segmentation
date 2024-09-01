function [Aa,bout,bb,ul] = mapacp_getborders(seg,bwdwm)

% example
% seg = MRIread2('atlas2mtr.nii'); seg = seg.vol;
% A = getweightedadjmatrix(seg);
%
% Alexandru V. Avram, Ph.D.
% NIH 2022

ul = unique(seg(seg>0)); % find all unique labels
bout = seg*0;
A = zeros(length(ul));
for k = 1:length(ul)
    % extract each label
    disp(k); tmp = seg*0; idx = find(seg==ul(k)); tmp(idx) = 1; 
    % find outside boundary layer (1 voxel thick)
    tmpp = bwperim(1-tmp); tmpp([1 end],:,:) = 0; tmpp(:,[1 end],:) = 0; tmpp(:,:,[1 end]) = 0; 
    % find bordering labels
    bvox = seg(find(tmpp)); ub = unique(bvox); ub(ub==0) = []; bb = seg.*tmpp;
    % remove bad connections
    for kb = 1:length(ub)
        if median(bwdwm(find(bb==ub(kb))))>7 | ub(kb)==0
            bb(find(bb==ub(kb))) = 0;
        end
    end
    % find new borders
    bbvox = seg(find(bb>0)); ubb = unique(bbvox);
    for kb = 1:length(ubb)
        nv = length(bbvox==ubb(kb)); % number of voxels
        A(k,find(ubb(kb)==ul)) = nv; % update adj mtx
    end
    bout = bout + (bb>0).*(seg>0);
end
bout = bout>0;

if ~issymmetric(A)
    Aa = A+A'; % make sure it's symmetric
else
    Aa = A;
end
%%
%for ii = 1:size(A,1)
%    for jj = 1:size(A,1)
%        Aa(ul(ii),ul(jj)) = A(ii,jj);
%    end
%end
