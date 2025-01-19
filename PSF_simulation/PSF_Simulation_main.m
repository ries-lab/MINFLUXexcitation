% This script simulates an excitation PSF for MINFLUX localization precision calculation. 
% Code written by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.
% This code utilizes a MATLAB software library, "Electromagnetic field in the focus region of a microscope objective", 
% written by Marcel Leutenegger. Request for the original package should be requested to marcel.leutenegger@alumni.epfl.ch.
% The library is included in this distribution.


%% Setting system parameters
% to automatically determine paths, this scirpt has to be run with the
% 'Run' functionality, not the 'Run section'. The sections below can be run
% with 'run section'.

clear
psimul=fileparts(mfilename('fullpath'));
addpath([psimul filesep 'library'])
if ~exist([psimul filesep 'simulated_data'],'dir')
    mkdir([psimul filesep 'simulated_data'])
end

sys = [];
out = [];    
opt.Et = 0; % Display Et during the calculation.
opt.Ef = 0;    % Display Ef during the calculation.
opt.calbar = 0;    % Display calculation progress bar.
opt.mem = 50000;
opt.pixSize = 2e-9; % Although the data size will be huge, pixel size of 2e-9 is recommended for fine calculation.
opt.radiusCanvas = 0.7e-6; % Lateral range.
opt.depthCanvas = 0.8e-6; % Axial range. 
opt.polAngle = 0.0*pi; % Linear polarization angle.
opt.phaseImage = false; % Show the phase image at the back focal plane.
opt.intImage = false;

% Setting specific parameters
[sys,out]=effInit_oil_exc(sys,out,opt);        % Assigning initial parameters. 
sys.NA=1.35;
sys.nm=1.406;
sys.ns=1.406;
sys.Mt=100;
sys.wa=7e-3;
sys.Na=200; 
sys.phaseImage = opt.phaseImage;
sys.intImage = opt.intImage;


%% Gaussian beam generation
disp('flat profile')
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.pl = 0; % Angle of the linear polarization.

sys.Zr=[2, 2, .1; 4,0,0.5];

sys.Ei = {'circular','zernike'};
% sys.Ei = {'circular'};
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);

figure(99); clf
hold all
plot((1:size(out.I, 3))*opt.pixSize, squeeze(out.I(1,round(y/2), :, round(z/2))))

simulation_flat_circular_xy = squeeze(out.I(1, :,:,round(z/2)));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_flat_circular_xy.mat'], 'simulation_flat_circular_xy')

maxgauss=simulation_flat_circular_xy(ceil(end/2),ceil(end/2));

sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.Ei = {'circular'};
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);

simulation_circular_xyz = squeeze(out.I(1, :,:,:));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_circular_xyz_pxlSize5nm.mat'], 'simulation_circular_xyz')

%% Vortex 2D donut.
disp('Vortex 2D donut')
sys.wa=7e-3;
sys.rz = 0.1e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.

sys.Ei = { 'phaseramp',  'circular'};    
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);
yaxis = linspace(0, y-1, y)*opt.pixSize;
zaxis = linspace(0, z-1, z)*opt.pixSize;
figure(99);
hold all
plot(yaxis, squeeze(out.I(1,round(y/2), :, round(z/2))))
hold all

simulation_vortex_circular_xy = squeeze(out.I(1, :,:,round(z/2)));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_vortex_circular_xy.mat'], 'simulation_vortex_circular_xy')

%% Calculation of PSF with halfmoon (-25 to 25 nm, xyz, pixel size 5 nm).
disp('bilobed')
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.pl = 0; % Angle of the linear polarization.
shiftmat = [-18.1, 0, 18.1]; % Shift of 25 nm, thus L = 50 nm.
matall = [];

for k = length(shiftmat):-1:1
    sys.delshift = deg2rad(shiftmat(k)); 
    sys.Ei = {'halfmoon', 'linear'};   
    
    %%% Calculating excitation PSF.
    % sys
    out=effField(sys,out, opt);      
    out=effIntensity(sys,out);
    
    [a,x,y,z] = size(out.I);
    matall(:,:,:, k)=squeeze(out.I(:,:,:,:));
end

simulation_halfmoon_linear_neg25to25nm_xyz = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_halfmoon_linear_neg25to25nm_xyz_pxlSize5nm.mat'], 'simulation_halfmoon_linear_neg25to25nm_xyz')

%% Calculation of PSF with halfmoon phase delays (displacement 0 to 150 nm).
disp('bilobed different L')
sys.wa=7e-3;
sys.rz = 0.01e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.
shiftmat = [0, 3.6,	7.3, 10.9, 14.5, 18.1, 21.8, 29, 36.3, 54.4, 72.6, 90.7, 108.9]; % Phase delay in degrees for displacements of 0, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 125, and 150 nm. 
% shiftmat = [-18.1, 0, 18.1]; % Shift of 25 nm, thus L = 50 nm.
matall = [];

figure(100);clf
for k = 1:length(shiftmat)
    sys.delshift = deg2rad(shiftmat(k));   
    sys.Ei = {'halfmoon', 'linear'};
    
%%% Calculating excitation PSF.
    % sys
    out=effField(sys,out, opt);      
    out=effIntensity(sys,out);
    [a,x,y,z] = size(out.I);
    yaxis = linspace(0, y-1, y)*opt.pixSize;
    zaxis = linspace(0, z-1, z)*opt.pixSize;
    hold all
    plot(yaxis, squeeze(out.I(1,round(y/2), :, round(z/2))))
    hold all
    matall(:,:,k)=squeeze(out.I(:,:,:,round(z/2)));
end

simulation_halfmoon_linear_0to150nm_xy = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_halfmoon_linear_0to150nm_xy.mat'], 'simulation_halfmoon_linear_0to150nm_xy')

%% Calculation of PSF with tophat (xy, no phase delay)
disp('tophat')
sys.wa=7e-3;
sys.rz = 0.1e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.delshift = 0;
sys.pl = 0; % Angle of the linear polarization.
sys.Ei = {'pishift', 'circular'};  
  
%%% Calculating excitation PSF.
% sys
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);

figure(99);
hold all
plot((1:size(out.I, 3))*opt.pixSize, squeeze(out.I(1,round(y/2), :, round(z/2))))
title('Intensity plot of different PSFs')
xlabel('Position in x-axis (m)')
ylabel('Intensity x-axis')
legend('Gaussian', '2D donut (Vortex)', '3D donut (Tophat)', 'Location', 'northeast')

simulation_tophat_circular_xy = squeeze(out.I(:,:,:,round(z/2)));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_tophat_circular_xy'], 'simulation_tophat_circular_xy')

%% Calculation of PSF with tophat ( -75 to 75 nm)
disp('tophat with different L')
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.

shiftmat = [20.84, 0, -20.84]; 
matall = [];
figure(98);clf
for k = length(shiftmat):-1:1
    sys.delshift = deg2rad(shiftmat(k));
    sys.Ei = {'pishift', 'circular'}; 

%%% Calculating excitation PSF.
    % sys
    out=effField(sys,out, opt);      
    out=effIntensity(sys,out);
    [a,x,y,z] = size(out.I);
    zaxis = linspace(0, z-1, z)*opt.pixSize;
    figure(98)
    hold all
    plot(zaxis, squeeze(out.I(1,round(y/2), round(x/2), :)))
    hold all
    matall(:,:,k)=squeeze(out.I(:,round(y/2),:,:)); 
end
figure(98);
title('3D donut PSF intensity plot at different phase')
xlabel('Position in z-axis (m)')
ylabel('Intensity z-axis')

simulation_tophat_circular_neg75to75nm_xz = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_tophat_circular_neg75to75nm_xz.mat'], 'simulation_tophat_circular_neg75to75nm_xz')

%% Calculation of PSF with tophat (- 75 to 75nm, pixel size 5 nm)
disp('tophat 3D')
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.pl = 0; % Angle of the linear polarization.

shiftmat = [20.84, 0, -20.84];
matall = [];

for k = length(shiftmat):-1:1
    sys.delshift = deg2rad(shiftmat(k));
    sys.Ei = {'pishift', 'circular'};   
%%% Calculating excitation PSF.
    out=effField(sys,out, opt);      
    out=effIntensity(sys,out);
    [a,x,y,z] = size(out.I);
    matall(:,:,:, k)=squeeze(out.I(:,:,:,:)); 
end

simulation_tophat_circular_neg75to75nm_xyz = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_tophat_circular_neg75to75nm_xyz_pxlSize5nm.mat'], 'simulation_tophat_circular_neg75to75nm_xyz')

%% Interferometric bi-lobe beam from Wolff, Scheiderer et al., Science 2023
disp('interferometric PSF')
% Beam diameter and positions are referred to a thesis by Dr. Tobias Engelhardt, where beam diameter 2mm, distance between two beams 4mm, back aperture 5.6mm.
sys.rz = 0.1e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
shiftx = [0.0023];
diam = [0.0023];
piston = [161.7, 180, 198.3];
sys.pl = 0.5*pi; % Angle of the linear polarization.
matall = [];
figure(102);clf;

for ii = 1:length(shiftx) % Loop for different beam shift amount.
    for kk = 1:length(diam) % Loop for different beam diameter.
        for jj = 1:length(piston) % Loop for different phase difference between beam pairs.
        
        clearvars outSum
        sys.Ei = {'offset','gauss', 'linear'};
        sys.wa = diam(kk);
        sys.ox = -shiftx(ii);
        sys.oy = 0;
        
        % Calculating the first beam.
        out=effField(sys,out, opt);      
        out=effIntensity(sys,out);
        [a,x,y,z] = size(out.E);
        Eright = out.E;
        
        % Calculating the second beam
        sys.Ei = {'piston', 'offset','gauss', 'linear'};
        sys.pistonphi = deg2rad(piston(jj)); % Phase delay for the other beam path.
        sys.ox = -sys.ox;
        sys.oy = -sys.oy;
        sys.pl = 0.5*pi;
        out=effField(sys,out, opt);      
        Eleft = out.E;
        
        % Calculating the sum of the two beams.
        outSum.E = Eright + Eleft;
        outSum=effIntensity(sys,outSum);
        outSum.I = outSum.I./2./141;
        figure(102);
        hold all;plot(((1:size(outSum.I, 3))-round(y/2))*opt.pixSize, squeeze((outSum.I(1,:,round(x/2),round(z/2))))./1)
        matall(:,:,jj) = squeeze(outSum.I(:,:,:,round(z/2)));
        end
    end
end
figure(102);
title('Interferometric PSF x-axis intensity plot at different phase')
xlabel('Position in x-axis (m)')
ylabel('Intensity x-axis')

simulation_gauss_linear_iMINFLUX_xyphi_L50nm = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_gauss_linear_iMINFLUX_xyphi_L50nm.mat'], 'simulation_gauss_linear_iMINFLUX_xyphi_L50nm')


%%  Calculation of PSF with halfmoon at wrong polarizations.
disp('bilobed with wrong polarizations')
sys.wa=7e-3;
sys.rz = 0.01e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
opt.pixSize = out.dr;
sys.pl = 0; % Angle of the linear polarization.
shiftmat = [0];
poldeg = linspace(0, 45, 46); % Degree of linear polarization orientation
matall = [];
figure(111);clf
for k = length(poldeg):-1:1
    sys.delshift = deg2rad(shiftmat);
    sys.pl = deg2rad(poldeg(k));
    sys.Ei = {'halfmoon', 'linear'};
    out=effField(sys,out, opt);      
    out=effIntensity(sys,out);
    [a,x,y,z] = size(out.I);
    xaxis = linspace(0, x-1, x)*opt.pixSize;
    hold all
    plot(xaxis, squeeze(out.I(1,round(y/2), :, round(z/2))))
    hold all
    matall(:,:,k)=squeeze(out.I(:,:,:,round(z/2)));
end
figure(111);
title('Halfmoon PSF x-axis intensity plot at wrong polarization angle')
xlabel('Polarization angle value from optimal')
ylabel('Intensity x-axis')

figure(110);clf
for kk = 1:size(matall, 3)
    hold all
    plot(kk-1, min(matall(round(y/2), round(x/2)-10:round(x/2)+10, kk))/maxgauss, '.')
end
title('Halfmoon PSF wrong polarization peak-to-minima')
xlabel('Polarization angle value from optimal')
ylabel('Zero-center intensity ratio')

simulation_halfmoon_linear_Pol_0to45deg_xyphi = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_halfmoon_linear_wrongPol_0to45deg_xyphi.mat'], 'simulation_halfmoon_linear_Pol_0to45deg_xyphi')
% end