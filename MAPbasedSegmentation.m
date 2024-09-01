%% Cortical Cytoarchitectonic Clustering
clear;clc;cd('D:\MacaqueCortParc')
x = load_untouch_nii('Res\atlas2mtr2024.nii'); a = single(x.img);

%% load data
load Res/MAPparams xpa xng xrtap xrtpp xfa xad xrd pial msp
%data = [xpa(:) xng(:) xrtap(:) xrtpp(:) xfa(:) xad(:) xrd(:)]; % MAP/DTI
data = [xpa(:) xng(:) xrtap(:) xrtpp(:)]; % MAP parameters
lidx = find(pial==1 & msp==1); ridx = find(pial==1 & msp==0);
[lhdata, lhc, lhs] = normalize(data(lidx,:)); [rhdata, rhc, rhs] = normalize(data(ridx,:));

%% LGMM clustering
% left hemisphere
nc = 14; options = statset('MaxIter',1000,'Display','iter'); 
lgmm = fitgmdist(lhdata,nc,'Options',options,'RegularizationValue',1e-4); 
lclust = cluster(lgmm,lhdata); lout = msp*0; lout(lidx) = lclust;

% right hemi segmentation
S.ComponentProportion = lgmm.ComponentProportion; S.mu = lgmm.mu; S.Sigma = lgmm.Sigma;
rgmm = fitgmdist(rhdata,nc,'Options',options,'Start',S,'RegularizationValue',1e-4); 
rclust = cluster(rgmm,rhdata); rout = msp*0; rout(ridx) = rclust; 

% morphological processing
llout = breakup_seg(lout,100); rrout = breakup_seg(rout,100);
lllout = merge_seg(llout,1e4); rrrout = merge_seg(rrout,1e4);

% matching with d99 layers
x = load_untouch_nii('Res\atlas2mtr2024layers.nii'); al = single(x.img);
x = load_untouch_nii('LAYNII\Final2024\rim2024_metric_equivol.nii'); rm = single(x.img);
x = load_untouch_nii('LAYNII\Final2024\rim2024_thickness.nii'); ct = single(x.img);
olt = match_seg(lllout,al.*msp,rm.*ct); ort = match_seg(rrrout,al.*(1-msp),rm.*ct);

% match clusters left-right hemis
xx = ort+olt; xort = ort; xolt = olt;
par = al.*(1-msp); pal = al.*msp;
T = readITKtable('Res/MAPACP_D99layers_crfs.txt'); Tval = table2array(T(:,2));
ort(1:length(Tval)) = Tval; olt(1:length(Tval)) = Tval;
par(1:length(Tval)) = Tval(end:-1:1); pal(1:length(Tval)) = Tval(end:-1:1);
[mr, chir, ptr, labr] = crosstab(ort(:),par(:)); mr(1,:) = []; mr(:,1) = []; labr(1,:) = []; 
[ml, chil, ptl, labl] = crosstab(pal(:),olt(:)); ml(1,:) = []; ml(:,1) = []; labl(1,:) = []; 
[rmm,rur,ruc] = matchpairs(mr*ml,5,'max');
for k = 1:length(rmm)
    xx(find(xort==rmm(k,1)))=rmm(k,2);
end
vol3d(xx) % saved as GMM2024Final.nii
