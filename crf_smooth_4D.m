function vv = crf_smooth_4D(vin,outx,ksz,sigsc)

% same as 3D but for 4D volume, e.g., dwis
%
% Alexandru V. Avram, Ph.D.
% NIH 2021
%

if nargin<4
    sigsc = [1/2 1/2 5];
end

if nargin<3
    ksz = 2;
end

sz = size(vin);
parfor k = 1:sz(4)
    disp(k); vv(:,:,:,k) = crf_smooth_3D(vin(:,:,:,k),outx,ksz,sigsc);
end