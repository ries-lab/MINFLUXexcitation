% CRB calculations for MINFLUX
% Jonas Ries, University of Vienna
% This follows https://github.com/stefani-lab/sml-ssi by Luciano Masullo
% Masullo, L.A., L.F. Lopez, and F.D. Stefani. 2022. A common framework for single-molecule localization using sequential structured illumination. Biophysical Reports. 2:100036. doi:10.1016/j.bpr.2021.100036.

function [sigma_CRB,sigma_CRBx,sigma_CRBy]=CRBMinflux(psfim, dx, dy,N,sbr,prior,s,sizeCRB)
    % Cramer-Rao Bound for a given SSI-SML experiment 
    % Input
    % ----------
    % psf: (K, size, size) array, experimental or simulated excitation beams
    % dx, dy: vector with shifts of the PSF in case only a single PSF image
    % is passed on
    % sbr : float, signal to background ratio
    % N : total number of photons
    % prior: type of a priori information, set to None if no prior is included
    % 'rough loc' is a model for a previous rough localization (gaussian-like)
    % s: σ parameter of the gaussian information prior, s=50 nm as default
    % rough localization, s >> L would be equivalent to no prior, s << L might
    % not be a realistic assumption
    % sizeCRB: size of the region for which the CRB is calculated
    % Output: σ_CRB,σ_CRBx,σ_CRBy, 
    
c=1;    % c: multiphoton order (e.g c=1 1-photon, c=2 2-photon)
sim=size(psfim);
if isempty(dx) || isempty(dy)% already PSFs at different positions passed on
    psf=psfim;
    K=sim(1);
else
    K=length(dy);
    dy=round(dy);dx=round(dx); %only integer shifts for now
    alld=[dy dx];
    maxshift=ceil(max(abs(alld(:))));
    clear psf
    psfim=squeeze(psfim); %allow for 2D image to be passed on
    for k=K:-1:1
        psf(k,:,:)=psfim(maxshift-dy(k)+1:end-maxshift-dy(k),maxshift-dx(k)+1:end-maxshift-dx(k));
    end
end
if nargin<8 || isempty(sizeCRB)
    sizeCRB=floor(size(psf,2:3)/2);
end

%size of the σ_CRB matrix in px and dimension d=2
size_matrix=size(psf); size_matrix(1)=[]; 
d = 2;
    
%XXX calculate in pixels, later convert to nm...
midp=ceil(size_matrix/2);
rangex=midp(1)-sizeCRB:midp(1)+sizeCRB; %for speed only calculate CRB within plotrange
rangey=midp(2)-sizeCRB:midp(2)+sizeCRB;

%initialize different arrays needed to compute σ_CRB, Σ_CRB and Fr
sigma_CRB = zeros(size_matrix(1), size_matrix(2));sigma_CRBx=sigma_CRB;sigma_CRBy=sigma_CRB;sigma_CRBam=sigma_CRB;
p=zeros(K,size_matrix(1), size_matrix(2)); dpdx=p;dpdy=p;
Fr_aux = zeros(K, d, d);
sbr_rel = sum(psf(:,rangex,rangey), 1);
sbr_rel=sbr_rel/sbr_rel(1,sizeCRB+1,sizeCRB+1); %normalize by central pixel?

sbr = sbr_rel * sbr;

%     # non-linearity due to multiphoton process of order c (c=1, linear regime)
psf = psf.^c;
%     # normalization term, sum of all psf intensities
norm_psf = sum(psf(:,rangex,rangey), 1);
for l=1:K
    p(l, rangex, rangey) = sbr./(sbr+1) .* psf(l,rangex,rangey)./norm_psf + (1./(sbr+1)) * (1/K);
    [dpdy(l, rangex, rangey), dpdx(l, rangex, rangey)] = gradient(squeeze(p(l, rangex, rangey)));
end

%     # compute relevant information for every (i, j) position within
%     plotrange
for i=midp(1)-sizeCRB:midp(1)+sizeCRB
    for j=midp(2)-sizeCRB:midp(2)+sizeCRB
        Fr_aux = zeros(d, d);
        for k=1:K
            A = [dpdx(k, i, j)^2, dpdx(k, i, j)*dpdy(k, i, j) ; dpdx(k, i, j)*dpdy(k, i, j), dpdy(k, i, j)^2];
            Fr_aux=Fr_aux+(1/p(k, i, j)) * A;
        end
        Fr=N*Fr_aux;
        if strcmp(prior,'rough loc')
             Fr=Fr+ diag([1/s^2, 1/s^2]);
        end
        SIGMA_CRB = inv(Fr);
        sigma_CRB(i, j) = sqrt((1/d) * (SIGMA_CRB(1,1)+SIGMA_CRB(2,2)));
        sigma_CRBx(i,j) = sqrt(SIGMA_CRB(1,1));
        sigma_CRBy(i,j) = sqrt(SIGMA_CRB(2,2));
    end
end
end