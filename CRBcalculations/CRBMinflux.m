function [sigma_CRB,sigma_CRBx,sigma_CRBy]=CRBMinflux(psfim, dx, dy,N,sbr,prior,s,sizeCRB)

c=1;
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
    % Cramer-Rao Bound for a given SSI-SML experiment 
    % Input
    % ----------
    % K : int, number of excitation beams
    % psf: (K, size, size) array, experimental or simulated excitation beams
    % SBR : float, signal to background ratio
    % px_nm : pixel of the grid in nm
    % size_nm : size of the grid in nm
    % N : total number of photons
    % c: multiphoton order (e.g c=1 1-photon, c=2 2-photon)
    % prior: type of a priori information, set to None if no prior is included
    % 'rough loc' is a model for a previous rough localization (gaussian-like)
    % s: σ parameter of the gaussian information prior, s=50 nm as default
    % rough localization, s >> L would be equivalent to no prior, s << L might
    % not be a realistic assumption
    % Output: σ_CRB, Σ_CRB, Fr, sbr_rel
    
%size of the σ_CRB matrix in px and dimension d=2
% size_matrix = int(size_nm/px_nm);
size_matrix=size(psf); size_matrix(1)=[]; 
d = 2;
    

%XXX calculate in pixels, later convert to nm...

%size of the (x,y) grid
% dx = px_nm;
% dy = px_nm;
        midp=ceil(size_matrix/2);
rangex=midp(1)-sizeCRB:midp(1)+sizeCRB;
rangey=midp(2)-sizeCRB:midp(2)+sizeCRB;

%initialize different arrays needed to compute σ_CRB, Σ_CRB and Fr
sigma_CRB = zeros(size_matrix(1), size_matrix(2));sigma_CRBx=sigma_CRB;sigma_CRBy=sigma_CRB;sigma_CRBam=sigma_CRB;
p=zeros(K,size_matrix(1), size_matrix(2)); dpdx=p;dpdy=p;
% Fr=zeros(d,d,size_matrix(1), size_matrix(2)); 
% Fr=zeros(d,d); 
% SIGMA_CRB = Fr;
% Fr_aux = zeros(K, d, d, size_matrix(1), size_matrix(2));
Fr_aux = zeros(K, d, d);

% sbr_rel = sum(psf, 1);
sbr_rel = sum(psf(:,rangex,rangey), 1);

        % sbr_rel = np.sum(psf, axis=0)/np.sum(psf, axis=0)[int(size/2), int(size/2)]
% sbr_rel=sbr_rel/sbr_rel(1,ceil(size_matrix(1)/2),ceil(size_matrix(2)/2)); %normalize by central pixel?
sbr_rel=sbr_rel/sbr_rel(1,sizeCRB+1,sizeCRB+1); %normalize by central pixel?

%     # SBR is computed as SBR(x, y), i.e. proportional to rel_sbr(x, y) instead of a constant
sbr = sbr_rel * sbr;

%     # non-linearity due to multiphoton process of order c (c=1, linear regime)
psf = psf.^c;





%     # normalization term, sum of all psf intensities
% norm_psf = sum(psf, 1);
norm_psf = sum(psf(:,rangex,rangey), 1);

for l=1:K
    % p(l, :, :) = sbr./(sbr+1) .* psf(l,:,:)./norm_psf + (1./(sbr+1)) * (1/K);
    p(l, rangex, rangey) = sbr./(sbr+1) .* psf(l,rangex,rangey)./norm_psf + (1./(sbr+1)) * (1/K);
%         # partial derivatives in x and y (minus sign for cartesian convention)
%         dpdy[i, :, :], dpdx[i, :, :] = np.gradient(p[i, :, :], -dy, dx)
    [dpdy(l, rangex, rangey), dpdx(l, rangex, rangey)] = gradient(squeeze(p(l, rangex, rangey)));
end
%
%     # compute relevant information for every (i, j) position
% for i=1:size_matrix
%     for j=1:size_matrix

for i=midp(1)-sizeCRB:midp(1)+sizeCRB
    for j=midp(2)-sizeCRB:midp(2)+sizeCRB
        Fr_aux = zeros(d, d);
        for k=1:K
 
            % A = [dpdx(k, i, j)^2, dpdx(k, i, j)*dpdy(k, i, j) ; dpdx(k, i, j)*dpdy(k, i, j), dpdy(k, i, j)^2];
            % Fr_aux(k, :, :, i, j) = (1/p(k, i, j)) * A;
            A = [dpdx(k, i, j)^2, dpdx(k, i, j)*dpdy(k, i, j) ; dpdx(k, i, j)*dpdy(k, i, j), dpdy(k, i, j)^2];
            % Fr_aux(k, :, :) = (1/p(k, i, j)) * A;
            Fr_aux=Fr_aux+(1/p(k, i, j)) * A;
        end
        % Fr(:,:) = N * (sum(Fr_aux, 1));
        Fr=N*Fr_aux;
        if strcmp(prior,'rough loc')
             Fr=Fr+ diag([1/s^2, 1/s^2]);
        end
        
        SIGMA_CRB = inv(Fr);
        % sigma_CRB(i, j) = sqrt((1/d) * trace(SIGMA_CRB));
        sigma_CRB(i, j) = sqrt((1/d) * (SIGMA_CRB(1,1)+SIGMA_CRB(2,2)));
        sigma_CRBx(i,j) = sqrt(SIGMA_CRB(1,1));
        sigma_CRBy(i,j) = sqrt(SIGMA_CRB(2,2));
        % eps=1e-12;
        % if det(Fr(:, :, i, j))<eps
        %     sigma_CRBx(i,j)=1/Fr(1,1,i,j);
        %     sigma_CRBy(i,j)=1/Fr(2,2,i,j);
        % end



        %weithed arithmetic mean:
        %<s^2>=1/sum(s_i^-2)

        % Ax=dpdx(k, i, j)^2;
        % Fr_aux1(k) = (1/p(k, i, j)) * A;
        % Fr(:, :, i, j) = N * sum(Fr_aux(:, :, :, i, j), 1);
        % SIGMA_CRB(:, :, i, j) = inv(Fr(:, :, i, j));
        % 
        % sigma_CRB(i,j) = sqrt(SIGMA_CRB(1,1,i,j));
        % 
        
    end
end
% sigma_CRBam = sqrt(1./sigma_CRBy.^2+sigma_CRBy.^2);


        %weithed arithmetic mean:
        %<s^2>=1/sum(s_i^-2)

end