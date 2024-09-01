function ort = match_seg(pl,al,bwd)

% match segmentation labels from a hemisphere to canonical layers 
% obtained from the symmetrical D99 atlas.
%
% Alexandru V. Avram, Ph.D.
% NIH 2022
%

% inputs
% pl - map segmentation (left/right hemi)
% al - symmetric d99_layers segmentation (left/right hemi)
% bwd - distance from white matter .e.g., rim_metric * rim_thickness

Nmx = 1e5;
pl(pl>0) = pl(pl>0)+Nmx;
[ml, chil, ptl, labl] = crosstab(al(:),pl(:)); ml(1,:) = []; ml(:,1) = []; labl(1,:) = []; 

% sort labels and recolor for left hemi
m = ml; lab = labl; th = [1 0.5 0.1 0.01]*1e2;
disp(['Filtering with threshold ',num2str(th(1))])
[MM,uR,uC] = matchpairs(m,th(1),'max');
for k = 1:length(MM)
    pl(find(pl==str2num(lab{MM(k,2),2}))) = str2num(lab{MM(k,1),1});
end
disp(['Assigned labels :',num2str(length(MM))])
disp(['Unassigned labels :',num2str(length(unique(pl(pl>Nmx))))])
mpl = m([MM(:,1);uR],[MM(:,2);uC]);

% set labelled regions to 0 and repeat
ppl = pl; ppl(ppl<Nmx) = 0; 
[ml, chil, ptl, labl] = crosstab(al(:),ppl(:)); ml(1,:) = []; ml(:,1) = []; labl(1,:) = []; 
m = ml; lab = labl;
disp(['Filtering with threshold ',num2str(th(2))])
[MM,uR,uC] = matchpairs(m,th(2),'max');
for k = 1:length(MM)
    ppl(find(ppl==str2num(lab{MM(k,2),2}))) = str2num(lab{MM(k,1),1});
end
disp(['Assigned labels :',num2str(length(MM))])
disp(['Unassigned labels :',num2str(length(unique(ppl(ppl>Nmx))))])
mppl = m([MM(:,1);uR],[MM(:,2);uC]);

% repeat
pppl = ppl; pppl(pppl<Nmx) = 0;
[ml, chil, ptl, labl] = crosstab(al(:),pppl(:)); ml(1,:) = []; ml(:,1) = []; labl(1,:) = []; 
m = ml; lab = labl;
disp(['Filtering with threshold ',num2str(th(3))])
[MM,uR,uC] = matchpairs(m,th(3),'max');
for k = 1:length(MM)
    pppl(find(pppl==str2num(lab{MM(k,2),2}))) = str2num(lab{MM(k,1),1});
end
disp(['Assigned labels :',num2str(length(MM))])
disp(['Unassigned labels :',num2str(length(unique(pppl(pppl>Nmx))))])
mpppl = m([MM(:,1);uR],[MM(:,2);uC]);

% last iteration
ppppl = pppl; ppppl(ppppl<Nmx) = 0;
[ml, chil, ptl, labl] = crosstab(al(:),ppppl(:)); ml(1,:) = []; ml(:,1) = []; labl(1,:) = []; 
m = ml; lab = labl;
disp(['Filtering with threshold ',num2str(th(4))])
[MM,uR,uC] = matchpairs(m,th(4),'max');
for k = 1:length(MM)
    ppppl(find(ppppl==str2num(lab{MM(k,2),2}))) = str2num(lab{MM(k,1),1});
end
disp(['Assigned labels :',num2str(length(MM))])
disp(['Unassigned labels :',num2str(length(unique(ppppl(ppppl>Nmx))))])

% fuse remaining labels with adjacent region sharing the largest boundary
rr = 0;
tmp = pl; tmp(tmp>Nmx) = 0; rr = tmp;
tmp = ppl; tmp(tmp>Nmx) = 0; rr = rr+tmp;
tmp = pppl; tmp(tmp>Nmx) = 0; rr = rr+tmp;
tmp = ppppl; tmp(tmp>Nmx) = 0; rr = rr+tmp;
tt = ppppl; tt(tt<Nmx) = 0; 
% finally assign last segments based on largest shared boundary
[Aa,bout,bb,ul] = mapacp_getborders(rr+tt,bwd);
rt = rr+tt; ort = rt;
idxl = unique(tt(tt>0));
rtul = unique(rt(rt>0));
for k = 1:length(idxl)
    currl = find(ul==idxl(k)); 
    tmp = Aa(currl,:); 
    if sum(tmp)~=0 
        [stmp, sx] = sort(tmp,'descend');
        if ismember(sx(1),find(rtul<Nmx))
            disp(['Assigning label ', num2str(idxl(k)), ' to value ', num2str(rtul(sx(1)))]);
            ort(find(ort==idxl(k))) = rtul(sx(1));
        end
    end
end
% assign to unlabeled for remaining voxels
ort(ort>Nmx) = 1;