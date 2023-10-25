% CRB calculations for MINFLUX
% Jonas Ries, University of Vienna
% This follows https://github.com/stefani-lab/sml-ssi by Luciano Masullo
% Masullo, L.A., L.F. Lopez, and F.D. Stefani. 2022. A common framework for single-molecule localization using sequential structured illumination. Biophysical Reports. 2:100036. doi:10.1016/j.bpr.2021.100036.
function [sigma_CRB,sigma_CRBx,sigma_CRBy,sigma_CRBz]=CRBMinflux3D(psf,sizeCRBxy)
% psf (K,size_x, size_y, size_z) Array describing PSF for K different
% probing positions
% sizeCRBxy: size of the region for which CRB is calculated
% output: sigma_CRB average, x, y, z normalized by sqrt(N), i.e., for a
% single photon
sim=size(psf);
K=sim(1);
if nargin<2 || isempty(sizeCRBxy)
    sizeCRBxy=floor(size(psf,2:3)/2);
end
size_matrix=size(psf); size_matrix(1)=[]; 
d = 3;

N=100; %normalize away later
sN=sqrt(N);
sizez=2; %defines slice size in z for gradient calculation
midp=ceil(size_matrix/2);
rangex=midp(1)-sizeCRBxy:midp(1)+sizeCRBxy;
rangey=midp(2)-sizeCRBxy:midp(2)+sizeCRBxy;
rangez=midp(3)-sizez:midp(3)+sizez;
indz=sizez+1;

%initialize different arrays needed to compute σ_CRB, Σ_CRB and Fr
sigma_CRB = zeros(length(rangex), length(rangey));sigma_CRBx=sigma_CRB;sigma_CRBy=sigma_CRB;sigma_CRBz=sigma_CRB;
p=zeros(K,length(rangex), length(rangey),length(rangez)); dpdx=p;dpdy=p;dpdz=p;

norm_psf = sum(psf(:,rangex,rangey,rangez), 1);

%calculate gradient
for l=1:K
    p(l, :, :,:) = psf(l,rangex,rangey,rangez)./norm_psf;
    [dpdy(l, :, :,:), dpdx(l, :,:,:), dpdz(l,:,:,:)] = gradient(squeeze(p(l, :, :,:)));
end

%     # compute relevant information for every (i, j) position at z=0
%     midplane
for i=1:length(rangex)
    for j=1:length(rangey)
        Fr_aux = zeros(d, d);
        for k=1:K
            A = [dpdx(k, i, j,indz)^2, dpdx(k, i, j,indz)*dpdy(k, i, j,indz), dpdx(k, i, j,indz)*dpdz(k, i, j,indz) ;...
                dpdx(k, i, j,indz)*dpdy(k, i, j,indz), dpdy(k, i, j,indz)^2, dpdy(k, i, j,indz)*dpdz(k, i, j,indz) ;...
                dpdx(k, i, j,indz)*dpdz(k, i, j,indz), dpdy(k, i, j,indz)*dpdz(k, i, j,indz),dpdz(k, i, j,indz)^2];
            Fr_aux=Fr_aux+(1/p(k, i, j,indz)) * A;
        end

        Fr=N*Fr_aux;
        SIGMA_CRB = inv(Fr);
        sigma_CRB(i,j) = sqrt((1/d) * (SIGMA_CRB(1,1)+SIGMA_CRB(2,2)+SIGMA_CRB(3,3)))*sN;
        sigma_CRBx(i,j) = sqrt(SIGMA_CRB(1,1))*sN;
        sigma_CRBy(i,j) = sqrt(SIGMA_CRB(2,2))*sN;
        sigma_CRBz(i,j) = sqrt(SIGMA_CRB(3,3))*sN;        
    end
end
end