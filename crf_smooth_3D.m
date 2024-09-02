function vv = crf_smooth_3D(vin,outx,ksz,sigsc)

% smooth volume based on crf
% inputs:   vin - volume (scalar)
%           outx - volume x 9 (reference frame) usually cortical
%           ksz - 1/2 size of kernel width (in voxels)
%           sigsc - standard deviations
% output:   vout - smoothed vin
%
% Alexandru V. Avram, Ph.D.
% NIH 2021
%

if nargin<4
    sigsc = [1/2 1/2 10];
end

if nargin<3
    ksz = 2;
end

vout = vin*0; sz = size(vin); isigs = ksz*sigsc;
[x,y,z] = meshgrid(-ksz:ksz,-ksz:ksz,-ksz:ksz); pts = [x(:) y(:) z(:)];
vinx = padarray(vin,[ksz ksz ksz]);
M = size(x,1)-1; N = size(y,1)-1; P = size(z,1)-1;
%Convolution
for kx = 1:sz(1)-M
    for ky =1:sz(2)-N
        for kz = 1:sz(3)-P
            if sum(outx(kx,ky,kz,:).^2)>0
                rpts = pts*reshape(squeeze(outx(kx,ky,kz,:)),[3 3]); 
                kernel = reshape(exp(-sum(rpts.*rpts.*repmat(isigs,[prod(size(x)) 1]),2)./(2*pi).^(3/2)/sqrt(1./prod(isigs))),size(x));
                tmp = vinx(kx:kx+M,ky:ky+N,kz:kz+P).*kernel./mean(kernel(:));
                vout(kx,ky,kz) = mean(tmp(:));
            end
        end
    end
    %kx 
end
vv = vin; vv(find(vout>0)) = vout(find(vout>0));