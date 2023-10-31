%% This function simulates an excitation PSF for MINFLUX localization precision calculation. 
% Code written by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.
% This code utilizes a MATLAB software library, "Electromagnetic field in the focus region of a microscope objective", 
% written by Marcel Leutenegger. Request for the original package should be requested to marcel.leutenegger@alumni.epfl.ch.
% The library is included in this distribution.

% function [out, sys, opt] = PSF_Simulation_main
%% Setting system parameters
clear
psimul=fileparts(mfilename('fullpath'));
addpath([psimul filesep 'library'])

sys = [];
out = [];    
opt.Et = 0; % Display Et during the calculation.
opt.Ef = 0;    % Display Ef during the calculation.
opt.calbar = 0;    % Display calculation progress bar.
opt.mem = 50000;
opt.pixSize = 5e-9; % Although the data size will be huge, pixel size of 2e-9 is recommended for fine calculation.
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
sys.wa=7e-3;
sys.rz = 0.1e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.

sys.Ei = {'circular'};
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);

figure(99); clf
title('Y axis profile')
hold all
plot((1:size(out.I, 3))*opt.pixSize, squeeze(out.I(1,round(y/2), :, round(z/2))))

simulation_silObj_planner_flat_circular_xy = squeeze(out.I(1, :,:,round(z/2)));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_flat_circular_xy.mat'], 'simulation_silObj_planner_flat_circular_xy')


sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.Ei = {'circular'}
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);

simulation_silObj_planner_circular_xyz = squeeze(out.I(1, :,:,:));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_circular_xyz_pxlSize5nm.mat'], 'simulation_silObj_planner_circular_xyz')



%% Vortex 2D donut.
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
figure(99);clf
hold all
plot(yaxis, squeeze(out.I(1,round(y/2), :, round(z/2))))
hold all

simulation_silObj_planner_vortex_circular_xy = squeeze(out.I(1, :,:,round(z/2)));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_vortex_circular_xy.mat'], 'simulation_silObj_planner_vortex_circular_xy')


%% Calculation of PSF with halfmoon (-25 to 25 nm, xyz, pixel size 5 nm).
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.pl = 0; % Angle of the linear polarization.
shiftmat = [-18.1, 0, 18.1]; % Shift of 25 nm, thus L = 50 nm.
matall = [];

for k = 1:length(shiftmat)

sys.delshift = deg2rad(shiftmat(k)); 
sys.Ei = {'halfmoon', 'linear'}   

%%% Calculating excitation PSF.
sys
out=effField(sys,out, opt);      
out=effIntensity(sys,out);

%%%% Plotting the 3D intensity
[a,x,y,z] = size(out.I);
matall(:,:,:, k)=squeeze(out.I(:,:,:,:));

end

simulation_silObj_planner_halfmoon_linear_neg25to25nm_xyz = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_halfmoon_linear_neg25to25nm_xyz_pxlSize5nm.mat'], 'simulation_silObj_planner_halfmoon_linear_neg25to25nm_xyz')


%% Calculation of PSF with halfmoon phase delays (displacement 0 to 150 nm).
sys.wa=7e-3;
sys.rz = 0.01e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.
shiftmat = [0, 3.6,	7.3, 10.9, 14.5, 18.1, 21.8, 29, 36.3, 54.4, 72.6, 90.7, 108.9]; % Phase delay in degrees for displacements of 0, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 125, and 150 nm. 
% shiftmat = [-18.1, 0, 18.1]; % Shift of 25 nm, thus L = 50 nm.
matall = [];

figure(99);clf
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

simulation_silObj_planner_halfmoon_linear_neg25to25nm_xy = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_halfmoon_linear_neg25to25nm_xy.mat'], 'simulation_silObj_planner_halfmoon_linear_neg25to25nm_xy')


%% Calculation of PSF with tophat (xy, no phase delay)
sys.wa=7e-3;
sys.rz = 0.01e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.
matall = [];
sys.Ei = {'pishift', 'circular'}  
  
%%% Calculating excitation PSF.
% sys
out=effField(sys,out, opt);      
out=effIntensity(sys,out);
[a,x,y,z] = size(out.I);


simulation_silObj_planner_tophat_circular_xy = squeeze(out.I(:,:,:,round(z/2)));
save([psimul, filesep, 'simulated_data', filesep, 'simulation_tophat_circular_xy'], 'simulation_silObj_planner_tophat_circular_xy')


%% Calculation of PSF with tophat ( -75 to 75 nm)
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
sys.pl = 0; % Angle of the linear polarization.

shiftmat = [20.84, 0, -20.84]; 
matall = [];
figure(98);clf
for k = 1:length(shiftmat)
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

simulation_silObj_planner_tophat_circular_neg75to75nm_xz = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_tophat_circular_neg75to75nm_xz.mat'], 'simulation_silObj_planner_tophat_circular_neg75to75nm_xz')


%% Calculation of PSF with tophat (- 75 to 75nm, pixel size 5 nm)
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
sys.pl = 0; % Angle of the linear polarization.

shiftmat = [20.84, 0, -20.84];
matall = [];

for k = 1:length(shiftmat)
    sys.delshift = deg2rad(shiftmat(k));
    sys.Ei = {'pishift', 'circular'}   

%%% Calculating excitation PSF.
    sys
    out=effField(sys,out, opt);      
    out=effIntensity(sys,out);
    [a,x,y,z] = size(out.I);
    matall(:,:,:, k)=squeeze(out.I(:,:,:,:)); 
end


simulation_silObj_planner_tophat_circular_neg75to75nm_xyz = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_tophat_circular_neg75to75nm_xyz_pxlSize5nm.mat'], 'simulation_silObj_planner_tophat_circular_neg75to75nm_xyz')




%% Interferometric bi-lobe beam from Wolff, Scheiderer et al., Science 2023
% Beam diameter and positions are referred to a thesis by Dr. Tobias Engelhardt, where beam diameter 2mm, distance between two beams 4mm, back aperture 5.6mm.
sys.rz = 0.1e-6;% Axial range.
out.dr = 2e-9;
out.dz = 2e-9;
shiftx = [0.0023];
diam = [0.0023];
piston = [161.7, 180, 198.3];
sys.pl = 0.5*pi; % Angle of the linear polarization.
matall = [];

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
        figure(102);clf;
        hold all;plot(((1:size(outSum.I, 3))-round(y/2))*opt.pixSize, squeeze((outSum.I(1,:,round(x/2),round(z/2))))./1)
        matall(:,:,jj) = squeeze(outSum.I(:,:,:,round(z/2)));
        end
    end
end
figure(102);
title('Interferometric PSF x-axis intensity plot at different phase')

simulation_gauss_linear_iMINFLUX_xyphi_L50nm = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_gauss_linear_iMINFLUX_xyphi_L50nm.mat'], 'simulation_gauss_linear_iMINFLUX_xyphi_L50nm')


%%  Calculation of PSF with halfmoon at wrong polarizations.
sys.wa=7e-3;
sys.rz = 0.01e-6;% Axial range.
out.dr = 5e-9;
out.dz = 5e-9;
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
    plot(kk-1, min(matall(round(y/2), round(x/2)-10:round(x/2)+10, kk))/max(matall(round(y/2), :, kk)), '.')
end
title('Halfmoon PSF wrong polarization peak-to-minima')
xlabel('Polarization angle value from optimal')
ylabel('Intensity ratio')

simulation_silObj_planner_halfmoon_linear_Pol_0to45deg_xyphi = matall;
save([psimul, filesep, 'simulated_data', filesep, 'simulation_halfmoon_linear_wrongPol_0to45deg_xyphi.mat'], 'simulation_silObj_planner_halfmoon_linear_Pol_0to45deg_xyphi')
% end