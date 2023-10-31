% CRB calculations for MINFLUX
% Jonas Ries, University of Vienna
% This follows https://github.com/stefani-lab/sml-ssi by Luciano Masullo
% Masullo, L.A., L.F. Lopez, and F.D. Stefani. 2022. A common framework for single-molecule localization using sequential structured illumination. Biophysical Reports. 2:100036. doi:10.1016/j.bpr.2021.100036.

clear 
pgit=fileparts(fileparts(fileparts(mfilename('fullpath'))));
smappath = [pgit filesep 'SMAP']; % this code uses helper function from SMAP, please specify the path explicitely if not founds
dirlist=genpath([smappath filesep 'shared']);addpath(dirlist)

simulationpath=[pgit filesep 'MINFLUXexcitation' filesep 'PSF_simulation' filesep 'simulated_data' filesep];

% important parameters
L=50; % nm, diameter (not radius) of scan pattern, only for PSFs which are not explicitely calculated for different positions
Lz=L*3;
bgoffsetrel=0.005; %background that is added to the PSF. It is defined in comparision to the maximum of an Airy PSF (flat phase, overfilled objective) with same the same integrated intensity
%for a vortex PSF at L=80 bgoffsetrel=0.01 corresponds to a center-frequncy-ratio (CFR) of 0.3

%additional parameters
shiftPSF=false; % no explicit calculation of PSFs (interference, bilobed, would work only for L=50 nm and Lz=150). Instead shift PSFs.
c=1; %multi-photon not considered here
prior='none'; % or 'rough loc'; then define
sprior=100; %pre-localization: width (sigma)  of prior
sbr=1000;% signal-to-background ratio. Assumed high, as we explicitely account for an imperfect contrast and a fluorescence background

%simulation data parameters
pixelsize=2; %nm, pixel size of simulations
pixelsize3D=5; %nm, coarser for 3D simulations
posval=[0, 5, 10, 15, 20, 25,30, 40, 50, 75, 100, 125, 150]; %simulations are saved at these positions
N=100; %photons, but dummy parameters, as this is canceled out by the normalization
indL=find(L/2==posval);
if isempty(L)
    error('no PSF calculated for given L size')
end

% plot parameters
plotrange=100; %nm;
maxcol=250; % color range for CRB plots in sigma * sqrt(N)
markerov='y.';
markerzoom='yo';

%% calculate CRB for defined L and various PSFs
%Airy PSF (flat phase, overfilled objective), used for calculating
%background offset for all PSFs
fileflat='simulation_flat_circular_xy.mat'; 
psfimflat=load([simulationpath fileflat]).simulation_flat_circular_xy;

gaussmaxint=max(psfimflat(:));
bgoffset=bgoffsetrel*gaussmaxint;

%Bilobed PSF for x and y localization
file='simulation_halfmoon_linear_0to150nm_xy'; % go back to all the values
psfimbl=load([simulationpath file]).simulation_halfmoon_linear_0to150nm_xy;
psfbl=makePSF2DL(psfimbl,indL); %Assembles PSF matrix from images
[sigma_CRB,sigma_CRBy,sigma_CRBx]=CRBMinflux(psfbl+bgoffset, [], [],N,sbr,prior,sprior,plotrange/pixelsize);

midpxy=floor((size(sigma_CRBx)+1)/2);
xall1=(-midpxy(1):midpxy(1))*pixelsize;
xall2=(-midpxy(2):midpxy(2))*pixelsize;

figure(100); clf
nx=4; ny=4; %subplots
subplot(nx,ny,1) %bilobed xy
imagesc(xall1,xall2,(psfimbl(:,:,1)+bgoffset)/gaussmaxint); colorbar
axis equal tight off
xlabel('x (nm)')
ylabel('y (nm)')
hold on
plot([midpxy(1)*pixelsize-110 midpxy(1)*pixelsize-10],[midpxy(2) midpxy(2)]*pixelsize-10,'w')
plotrect(gca,[-plotrange -plotrange plotrange plotrange],'w')
plot([-1 0 1]*L/2,[0 0 0]*L/2,markerov,'MarkerSize',3)
intenL=squeeze(psfimbl(midpxy(1),midpxy(2),indL))/gaussmaxint;% intensity at scan position
title(['bisected, I(L/2) = ' num2str(intenL,2)])
subplot(nx,ny,ny+1) %bilobed crb
hold off
plotPSF(sigma_CRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm)')
hold on
plot([-1 0 1 0 0 0]*L/2,[0 0 0 -1 0 1]*L/2,markerzoom,'MarkerSize',3)

%2D donut with a vortex phase
filevortex='simulation_vortex_circular_xy.mat';
psfimvortex=load([simulationpath filevortex]).simulation_vortex_circular_xy;
%scan pattern: triangle + center
patternx3(:,1)=[0,1,-0.5,-0.5];
patterny3(:,1)=[0,0,sind(120),-sind(120)];
[sigma_vCRB,sigma_vCRBy,sigma_vCRBx]=CRBMinflux(psfimvortex+bgoffset, patternx3*L/2/pixelsize, patterny3*L/2/pixelsize,N,sbr,prior,sprior,plotrange/pixelsize);
subplot(nx,ny,3)
imagesc(xall2,xall1,(psfimvortex+bgoffset)/gaussmaxint); colorbar
axis equal tight off
xlabel('x (nm)')
ylabel('z (nm)')
intenL=psfimvortex(midpxy(1)+round(L/pixelsize/2),midpxy(2),1)/gaussmaxint;
title(['vortex, I(L/2) = ' num2str(intenL,2)])
hold on
plot([midpxy(1)*pixelsize-110 midpxy(1)*pixelsize-10],[midpxy(2) midpxy(2)]*pixelsize-10,'w')
plotrect(gca,[-plotrange -plotrange plotrange plotrange],'w')
plot(patternx3*L/2,patterny3*L/2,markerov,'MarkerSize',3)
hold off
subplot(nx,ny,3+ny)
plotPSF(sigma_vCRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm)')
hold on
plot(patternx3*L/2,patterny3*L/2,markerzoom,'MarkerSize',3)

% 3D donut with tophat phase, 2D scanning
fileth='simulation_tophat_circular_xy.mat';
psfimth=load([simulationpath fileth]).simulation_tophat_circular_xy;
[sigma_thCRB,sigma_thCRBy,sigma_thCRBx]=CRBMinflux(psfimth+bgoffset, patternx3*L/2/pixelsize, patterny3*L/2/pixelsize,N,sbr,prior,sprior,plotrange/pixelsize);

subplot(nx,ny,9)
imagesc(xall2,xall1,(psfimth+bgoffset)/gaussmaxint); colorbar
axis equal tight off
xlabel('x (nm)')
ylabel('z (nm)')
intenL=psfimth(midpxy(1)+round(L/pixelsize/2),midpxy(2),ceil(end/2))/gaussmaxint;
title(['top hat, I(L/2) = ' num2str(intenL,2)])
hold on
plot([midpxy(1)*pixelsize-110 midpxy(1)*pixelsize-10],[midpxy(2) midpxy(2)]*pixelsize-10,'w')
plotrect(gca,[-plotrange -plotrange plotrange plotrange],'w')
plot(patternx3*L/2,patterny3*L/2,markerov,'MarkerSize',3)

subplot(nx,ny,ny+9);
plotPSF(sigma_thCRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm)')
hold on
plot(patternx3*L/2,patterny3*L/2,markerzoom,'MarkerSize',3)

% flat phase, Airy PSF with overfille objective
L2x2=300;%larger scanning 
pattern2x2x=[-1 -1 1 1]; %2 x 2 scan pattern
pattern2x2y=[-1 1 -1 1];
pattern3x3x=[-1 -1 -1 0 0 0 1 1 1]; %3 x 3 scan pattern
pattern3x3y=[-1 0 1 -1 0 1 -1 0 1];
patternGx=pattern3x3x; %choose which scan pattern to display
patternGy=pattern3x3y;
plotrangeG=250; %larger range for Airy PSF
[sigma_fCRB,sigma_fCRBy,sigma_fCRBx]=CRBMinflux(psfimflat+bgoffset, patternGx*L2x2/2/pixelsize, patternGy*L2x2/2/pixelsize,N,sbr,prior,sprior,plotrangeG/pixelsize);

subplot(nx,ny,11)
imagesc(xall2,xall1,(psfimflat+bgoffset)/gaussmaxint); colorbar
axis equal tight off
xlabel('x (nm)')
ylabel('y (nm)')
title(['scanning size L = ' num2str(L2x2) ' nm'])
hold on
plot([midpxy(1)*pixelsize-110 midpxy(1)*pixelsize-10],[midpxy(2) midpxy(2)]*pixelsize-10,'w')
plotrect(gca,[-plotrangeG -plotrangeG plotrangeG plotrangeG],'w')
plot(patternGx*L2x2/2,patternGy*L2x2/2,markerov,'MarkerSize',3)

subplot(nx,ny,ny+11)
plotPSF(sigma_fCRB*sqrt(N),plotrangeG,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm)')
hold on
plot(patternGx*L2x2/2,patternGy*L2x2/2,markerzoom,'MarkerSize',3)

% interferometric psf from Wolff et al
patternx=[0 0 0]; patterny=[-1 0 1]; %1D scan pattern
fileif='simulation_gauss_linear_iMINFLUX_xyphi_L50nm.mat';
psfimif=load([simulationpath fileif]).simulation_gauss_linear_iMINFLUX_xyphi_L50nm;
psfimif=psfimif/sqrt(2); %wrong normalization during simulation
if shiftPSF
    psfif=makePSF2DLmirror(psfimif(:,:,2),L/pixelsize);
    [sigma_ifCRB,sigma_ifCRBy,sigma_ifCRBx]=CRBMinflux(psfif+bgoffset, [], [],N,sbr,prior,sprior,plotrange/pixelsize);
else
    %make PSF for x and y scanning, calculated for a phase shift corresponding
    %to L=50. Virutally this was the same as shifting the PSF in space (above)
    clear psfifd
    psfifd(1,:,:)=psfimif(:,:,1);
    psfifd(2,:,:)=psfimif(:,:,2);
    psfifd(3,:,:)=psfimif(:,:,3);
    psfifd(4,:,:)=psfimif(:,:,1)';
    psfifd(5,:,:)=psfimif(:,:,2)';
    psfifd(6,:,:)=psfimif(:,:,3)';
    [sigma_ifCRB,sigma_ifCRBy,sigma_ifCRBx]=CRBMinflux(psfifd+bgoffset, [], [],N,sbr,prior,sprior,plotrange/pixelsize);
end
subplot(nx,ny,12)
imagesc(xall2,xall1,(psfimif(:,:,2)+bgoffset)/gaussmaxint); colorbar
axis equal tight off
xlabel('x (nm)')
ylabel('z (nm)')
intenL=psfimif(midpxy(1)+round(L/pixelsize/2),midpxy(2),2)/gaussmaxint;
title(['Interferometric, I(L/2) = ' num2str(intenL,2)])
hold on
plot([midpxy(1)*pixelsize-110 midpxy(1)*pixelsize-10],[midpxy(2) midpxy(2)]*pixelsize-10,'w')
plotrect(gca,[-plotrange -plotrange plotrange plotrange],'w')
plot(patternx*L/2,patterny*L/2,markerov,'MarkerSize',3)
plot(patterny*L/2,patternx*L/2,markerov,'MarkerSize',3)

subplot(nx,ny,ny+12)
hold off
plotPSF(sigma_ifCRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm)')
hold on
plot(patternx*L/2,patterny*L/2,markerzoom,'MarkerSize',3)
plot(patterny*L/2,patternx*L/2,markerzoom,'MarkerSize',3)

%% 3D CRBs
%flat for normalization of background:
filef='simulation_circular_xyz_pxlSize5nm.mat';
psfimf=load([simulationpath filef]).simulation_circular_xyz;
cslice=psfimf(:,:,ceil(end/2));
gaussmaxint3D=max(cslice(:));
bgoffset3D=bgoffsetrel*gaussmaxint3D;

fileth3D='simulation_tophat_circular_neg75to75nm_xyz_pxlSize5nm.mat';
psfimth3D=load([simulationpath fileth3D]).simulation_tophat_circular_neg75to75nm_xyz;

fileb3D='simulation_halfmoon_linear_neg25to25nm_xyz_pxlSize5nm.mat';
psfimb3D=load([simulationpath fileb3D]).simulation_halfmoon_linear_neg25to25nm_xyz;

% 3D scannign of tophat
dx=[-1 -1 1 1 0 0 0];
dy=[-1 1 -1 1 0 0 0];
dz=[0 0 0 0 -1 1 0];

% same scan pattern as for bisected PSF, reflects also the longer time
% spent in center in Abberior MINFLUX
dx=[1 0 -1 0 0 0  0 0 0];
dy=[ 0 0 0  1 0 -1 0 0 0];
dz=[ 0 0 0  0 0 0  1 0 -1];

%2D bilobed + z with 3D donut, shifted for any L
if shiftPSF
    psf3D=makeshiftedPSF3D(psfimb3D(:,:,:,2), [-1 0 1]/pixelsize3D*L/2, [0 0 0],[0 0 0]);
    psf3D(4:6,:,:,:)=permute(psf3D,[1,3,2,4]); %rotate by 90° for y
    psfz=makeshiftedPSF3D(psfimth3D(:,:,:,2), [0 0 0], [0 0 0], [-1 0 1]/pixelsize3D*Lz/2);
    psf3D(7:9,:,:,:)=psfz;
    
    
    [sigma_3DxyCRB,sigma_3DxyCRBx,sigma_3DxyCRBy,sigma_3DxyCRBz]=CRBMinflux3D(psf3D+bgoffset3D,plotrange/pixelsize3D); %x-y slice
    psf3Dxz=permute(psf3D,[1,2,4,3]); %x-z slice
    [sigma_3DxzCRB,sigma_3DxzCRBx,sigma_3DxzCRBz,sigma_3DxzCRBy]=CRBMinflux3D(psf3Dxz+bgoffset3D,plotrange/pixelsize3D);
else
    % use PSF with calculated phase shifts (L=50 nm in xy, L=150 nm in z)
    psf3Dc=permute(psfimb3D,[4,1,2,3]);
    psf3Dc(4:6,:,:,:)=permute(psfimb3D,[4,2,1,3]);
    psf3Dc(7:9,:,:,:)=permute(psfimth3D,[4,1,2,3]);
    
    [sigma_3DCRB,sigma_3DCRBx,sigma_3DCRBy,sigma_3DCRBz]=CRBMinflux3D(psf3Dc+bgoffset3D,plotrange/pixelsize3D);
    psf3Dcxz=permute(psf3Dc,[1,2,4,3]); %make x-z slice
    [sigma_3DxzCRB,sigma_3DxzCRBx,sigma_3DxzCRBz,sigma_3DxzCRBy]=CRBMinflux3D(psf3Dcxz+bgoffset3D,plotrange/pixelsize3D);
end

subplot(nx,ny,2)
xallz=(ceil(-size(psfimb3D,1)/2):floor(size(psfimb3D,1)/2))*pixelsize3D;
xallz2=(ceil(-size(psfimb3D,3)/2):floor(size(psfimb3D,3)/2))*pixelsize3D;
imagesc(xallz,xallz2,squeeze(psfimth3D(ceil(end/2),:,:,2)+bgoffset)'/gaussmaxint); colorbar
axis equal tight off
xlabel('x (nm)')
ylabel('z (nm)')
hold on
plot(dx*L/2,dz*Lz/2, markerov,'MarkerSize',3)
plotrect(gca,[-plotrange -plotrange plotrange plotrange],'w')
title('tophat')

subplot(nx,ny,2+ny)
plotPSF(sigma_3DxzCRB',plotrange,pixelsize3D,maxcol,'bisect + tophat 3D xz')
hold on
plot(dx*L/2,dz*Lz/2, markerzoom,'MarkerSize',3)


%3D donut, tophat profile, 3D scanning
psf3Dth=makeshiftedPSF3D(psfimth3D(:,:,:,2), dx/pixelsize3D*L/2, dy/pixelsize3D*L/2,dz/pixelsize3D*Lz/2);
% [sigma_thCRB,sigma_thCRBx,sigma_thCRBy,sigma_thCRBz]=CRBMinflux3D(psf3Dth+bgoffset3D,plotrange/pixelsize3D); %x-y slice
psf3Dxzth=permute(psf3Dth,[1,2,4,3]); %now x-z slice
[sigma_thxzCRB,sigma_thxzCRBx,sigma_thxzCRBz,sigma_thxzCRBy]=CRBMinflux3D(psf3Dxzth+bgoffset3D,plotrange/pixelsize3D);

subplot(nx,ny,10+ny)
plotPSF(sigma_thxzCRB',plotrange,pixelsize3D,maxcol,'tophat 3D xz')
hold on
plot(dx*L/2,dz*Lz/2, 'ko','MarkerSize',3)


%% intensity vs L plot
%load x-z image of tophat
filen3Dthxz='simulation_tophat_circular_neg75to75nm_xz.mat';
psfimz=load([simulationpath filen3Dthxz]).simulation_tophat_circular_neg75to75nm_xz;
midpz=ceil(size(psfimz)/2);

dL=max(posval)/pixelsize;
intenLbsxy(1,:)=squeeze(psfimbl(midpxy(1),midpxy(2):midpxy(2)+dL,1))/gaussmaxint;
intenLthz(1,:)=squeeze(psfimz(midpz(1),midpz(2):midpz(2)+dL,ceil(end/2)))/gaussmaxint;
intenLthxy(1,:)=squeeze(psfimth(midpxy(1):midpxy(1)+dL,midpxy(2)))/gaussmaxint;
intenLv(1,:)=squeeze(psfimvortex(midpz(1):midpz(1)+dL,midpxy(2)))/gaussmaxint;
intenLif(1,:)=squeeze(psfimif(midpz(1):midpz(1)+dL,midpxy(2),2))/gaussmaxint;

figure(101); clf
subplot(2,2,1)
xv=(0:pixelsize:max(posval))*2;
plot(xv,intenLbsxy)
hold on
plot(xv,intenLthz)
plot(xv,intenLthxy)
plot(xv,intenLv)
plot(xv,intenLif)
xlabel('L (nm)')
ylabel('I(L) norm')

ltxt={'bisect: x','tophat: z', 'tophat: x', 'vortex: x', 'interferometric: x'};
legend(ltxt)

subplot(2,2,2)
plot(xv(2:end),intenLbsxy(2:end)./intenLbsxy(2:end))
hold on
plot(xv(2:end),intenLbsxy(2:end)./intenLthz(2:end))
plot(xv(2:end),intenLbsxy(2:end)./intenLthxy(2:end))
plot(xv(2:end),intenLbsxy(2:end)./intenLv(2:end))
plot(xv(2:end),intenLbsxy(2:end)./intenLif(2:end))

xlabel('L (nm)')
ylabel('I(L)/I_b(L)')
title('Laser intensity compared to bisected')
legend(ltxt)

%% CRB at center for different L
subplot(2,2,3)
bgs=[0 0.005 0.01];
styles={'-','--','-.'};
plot(0,0,['k' styles{1}],0,0,['k' styles{2}],0,0,['k' styles{3}]);%,0,0,['k' styles{4}])
hold on
clear profb profth profv lt
for k=1:length(bgs)
    [sL2D,sLy2D,sLx2D,~,profb(:,:,k)]=getCRBL(psfimbl, N, sbr,bgs(k)*gaussmaxint, prior, sprior,plotrange/pixelsize,@makePSF2DL);
    % [sLz,sLyz,sLxz,psfz,prof2]=getCRBL(psfimz, N, sbr,bgs(k)*gaussmaxint, prior, sprior,plotrange/pixelsize,@makePSFz,length(posval));
    [sLth2D,sLyth2D,sLxth2D,~,profth(:,:,k)]=getCRBLshift(psfimth+bgs(k)*gaussmaxint, N, sbr, prior, sprior,patternx3, patterny3,posval/pixelsize,plotrange/pixelsize);
    [sLv2D,sLyv2D,sLxv2D,~,profv(:,:,k)]=getCRBLshift(psfimvortex+bgs(k)*gaussmaxint, N, sbr, prior, sprior,patternx3, patterny3,posval/pixelsize,plotrange/pixelsize);
    [sigma_fCRB,sigma_fCRBy,sigma_fCRBx]=CRBMinflux(psfimflat+bgs(k)*gaussmaxint, pattern2x2x*L2x2/2/pixelsize, pattern2x2y*L2x2/2/pixelsize,N,sbr,prior,sprior,plotrange/pixelsize);
    plot(posval(2:end)*2,sL2D(2:end)*sqrt(N)*pixelsize,['b' styles{k}])
    hold on
    plot(posval(2:end)*2,sLth2D(2:end)*sqrt(N)*pixelsize,['k' styles{k}])
    plot(posval(2:end)*2,sLv2D(2:end)*sqrt(N)*pixelsize,['r' styles{k}])
    plot(posval([2 end])*2,sigma_fCRB(ceil(end/2),ceil(end/2))*sqrt(N)*pixelsize*[1 1],['g' styles{k}])
    lt{k}=['bg: ' num2str(bgs(k))];
end

xlabel('L')
ylabel('$$\sigma \cdot \sqrt{N}$$ (nm)' ,'Interpreter','latex')
title('\sigma vs L')
lt(end+1:end+4)={'bisect xy 2D','tophat xy 2D','vortex 2D','Gauss L=300 nm'};
legend(lt,"Location","southeast")
xlim([0 200])
ylim([0 150])

%% field of view 
%2D maps along x
indLs=[5 8]+1; % which L to sample
posval(indLs)*2
nxv=(0:pixelsize:plotrange)';

subplot(2,2,4)
clear lt
styles={'-','--','-.',':'};
plot(0,0,['k' styles{1}],0,0,['k' styles{2}],0,0,['k' styles{3}]);%,0,0,['k' styles{4}])
hold on
for l=1:length(indLs)
    plot(nxv,profb(:,indLs(l),2)*sqrt(N)*pixelsize,[styles{l} 'b'])
    plot(nxv,profth(:,indLs(l),2)*sqrt(N)*pixelsize,[styles{l} 'k'])
    plot(nxv,profv(:,indLs(l),2)*sqrt(N)*pixelsize,[styles{l} 'r'])
    lt{l}=['L = ' num2str(posval(indLs(l))*2)];
end
lt(end+1:end+3)={'bisected','tophat','vortex'};
legend(lt)
xlabel('Position in FoV (nm)')
ylabel('$$\sigma_{2D} \cdot \sqrt{N}$$ (nm)','Interpreter','latex')
title('CRB vs FoV, bg 0.005')
ylim([0 200])


%% try more than 3 steps
%do for 2D at least
psfbs=psfbl(1,:,:);
pattern3=[-1 0 1];
pattern5=[-1 -0.5 0 0.5 1];
pattern21=-1:0.1:1;

psf3=makeshiftedPSF(psfbs, pattern3*L/2/pixelsize, pattern3*0);
psf3(end+1:2*end,:,:)=permute(psf3,[1 3 2]);
psf5=makeshiftedPSF(psfbs, pattern5*L/2/pixelsize, pattern5*0);
psf5(end+1:2*end,:,:)=permute(psf5,[1 3 2]);
psf21=makeshiftedPSF(psfbs, pattern21*L/2/pixelsize, pattern21*0);
psf21(end+1:2*end,:,:)=permute(psf21,[1 3 2]);

[sigma_3CRB,sigma_3CRBy,sigma_3CRBx]=CRBMinflux(psf3+bgoffset, [], [],N,sbr,prior,sprior,plotrange/pixelsize);
[sigma_5CRB,sigma_5CRBy,sigma_5CRBx]=CRBMinflux(psf5+bgoffset, [], [],N,sbr,prior,sprior,plotrange/pixelsize);
[sigma_21CRB,sigma_21CRBy,sigma_21CRBx]=CRBMinflux(psf21+bgoffset, [], [],N,sbr,prior,sprior,plotrange/pixelsize);

midps=ceil(size(sigma_3CRB)/2);

figure(103);clf
subplot(2,2,1)
plotPSF(sigma_3CRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm), 3 steps')
subplot(2,2,2)
plotPSF(sigma_5CRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm), 5 steps')
subplot(2,2,3)
plotPSF(sigma_21CRB*sqrt(N),plotrange,pixelsize,maxcol,'$$\sigma_{2D} \cdot \sqrt{N}$$ (nm), 21 steps')

subplot(2,2,4)
nxv=(1:midps(2))*pixelsize;
plot(nxv,sigma_3CRB(midps(1),ceil(end/2):end)*pixelsize*sqrt(N))
hold on
plot(nxv,sigma_5CRB(midps(1),ceil(end/2):end)*pixelsize*sqrt(N))
plot(nxv,sigma_21CRB(midps(1),ceil(end/2):end)*pixelsize*sqrt(N))
legend('3 steps','5 steps','21 steps')
xlim([0 100])


function plotPSF(imp,plotrangein,pixelsize,maxcol,titlet)
midp=floor((size(imp)+1)/2);
plotrange=round(plotrangein/pixelsize);
imp(isinf(imp))=maxcol;
imp(isnan(imp))=maxcol;
plotim=imp(midp(1)-plotrange:midp(1)+plotrange,midp(2)-plotrange:midp(2)+plotrange)*pixelsize;
xplot=(-plotrange:plotrange)*pixelsize;
imagesc(xplot,xplot,plotim,[0 maxcol]); colorbar
axis equal tight off
cols=parula;
cols(end,:)=[1,1,1];
colormap(gca,cols)
title([titlet ' , min= ' num2str(min(plotim(:)),3)],'Interpreter','latex')
xlabel('x (nm)')
ylabel('y (nm)')

hold on
xline=[plotrange-100/pixelsize plotrange]-1;
yline=(plotrange-1)*[1 1];
plot(xline*pixelsize,yline*pixelsize,'k')
hold off
end

function psf=makePSF2DLmirror(p1,Lpix)
p1=imtranslate(p1,[0.1,0.1],"cubic"); % avoid indeterministic vvalues at zero
psfim=p1;
psfim_flipped=p1';
dy=[-1 0 1]*Lpix/2;
dx=[0 0 0]*Lpix/2;
psf=makeshiftedPSF(psfim, dx, dy);
psf(end+1:end+3,:,:)=makeshiftedPSF(psfim_flipped, dy, dx);
end

function psf=makePSF2DL(p1,indL)
p1=imtranslate(p1,[0.1,0.1],"cubic"); % avoid indeterministic values at zero
psfim=p1;
if nargin>1
    psfim_mirrored=p1(:,end:-1:1,:);
    psf(1,:,:)=psfim(:,:,1);
    psf(2,:,:)=psfim(:,:,indL);
    psf(3,:,:)=psfim_mirrored(:,:,indL);
    psf(4,:,:)=psfim(:,:,1)';
    psf(5,:,:)=psfim(:,:,indL)';
    psf(6,:,:)=psfim_mirrored(:,:,indL)';
else
    psf(1,:,:)=psfim(:,:,1);
    psf(2,:,:)=psfim(:,:,2);
    psf(3,:,:)=psfim(:,:,3);
    psf(4,:,:)=psfim(:,:,1)';
    psf(5,:,:)=psfim(:,:,2)';
    psf(6,:,:)=psfim(:,:,3)';    
end
end

function psfz=makePSFz(pz,indL)
pz=imtranslate(pz,[0.1,0.1],"cubic"); % avoid indeterministic vvalues at zero
psfimz=pz;
midp=ceil(size(pz,3)/2);
psfz(1,:,:)=psfimz(:,:,midp);
psfz(2,:,:)=psfimz(:,:,midp+indL-1);
psfz(3,:,:)=psfimz(:,:,midp-indL+1);
end

function [sL,sLy,sLx,psf,prof2]=getCRBL(psfim, N, sbr,backgroundoffset, prior, sprior,plotrangepix,makePSFfunction,numL)
if nargin<9
    numL=size(psfim,3);
end
for l=2:numL
    psf=makePSFfunction(psfim,l)+backgroundoffset;
    [sigma_LCRB,sigma_LCRBy,Lsigma_LCRBx]=CRBMinflux(psf, [], [],N,sbr,prior,sprior,plotrangepix);
    midpxy=ceil(size(sigma_LCRB)/2);
    sL(l)=sigma_LCRB(midpxy(1),midpxy(2));sLy(l)=sigma_LCRBy(midpxy(1),midpxy(2));sLx(l)=Lsigma_LCRBx(midpxy(1),midpxy(2));
    prof2(:,l)=makecrbprofile(sigma_LCRB,plotrangepix);
end
end

function [sL,sLy,sLx,psfim,prof2]=getCRBLshift(psfim, N, sbr,prior, sprior,patternxth, patternyth,posvalpix,plotrangepix)
for l=2:length(posvalpix)
    [sigma_LCRB,sigma_LCRBy,Lsigma_LCRBx]=CRBMinflux(squeeze(psfim), patternxth*posvalpix(l), patternyth*posvalpix(l),N,sbr,prior,sprior,plotrangepix);
    midpxy=ceil(size(sigma_LCRB)/2);
    sL(l)=sigma_LCRB(midpxy(1),midpxy(2));sLy(l)=sigma_LCRBy(midpxy(1),midpxy(2));sLx(l)=Lsigma_LCRBx(midpxy(1),midpxy(2));
    prof2(:,l)=makecrbprofile(sigma_LCRB,plotrangepix);
end
end

function psf=makeshiftedPSF3D(psfim, dx, dy,dz)
   K=length(dy);
    dy=round(dy);dx=round(dx); dz=round(dz);%only integer shifts for now
    alld=[dy dx];
    maxshift=ceil(max(abs(alld(:))));
    psfim=squeeze(psfim); %allow for 3D image to be passed on
    psf=zeros(K,size(psfim,1),size(psfim,2),size(psfim,3));
    maxdz=max(abs(dz));
    for k=K:-1:1
        psf(k,maxshift+1:end-maxshift,maxshift+1:end-maxshift,maxdz+1:end-maxdz)=psfim(maxshift-dy(k)+1:end-maxshift-dy(k),maxshift-dx(k)+1:end-maxshift-dx(k),maxdz-dz(k)+1:end-maxdz-dz(k));
    end
end

function psf=makeshiftedPSF(psfim, dx, dy)
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

function imout=makecrbprofile(imin,wi)
imout(1,:)=imin(ceil(end/2),ceil(end/2):ceil(end/2)+wi);
imout(1,:)=imout+imin(ceil(end/2),ceil(end/2):-1:ceil(end/2)-wi);
imout(1,:)=imout+imin(ceil(end/2):ceil(end/2)+wi,ceil(end/2))';
imout(1,:)=imout+imin(ceil(end/2):-1:ceil(end/2)-wi,ceil(end/2))';
imout=imout/4;
end
