% Code modified by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.



%[sys,out]=effInit(sys,out)
%--------------------------
%
%Initialize the computation of the focus field.
%
% sys    System data
%  .Da      Aperture diameter [m]                     {6.5mm}
%  .da      Aperture sampling [m]                     {inf:auto}
%  .Ei      Field functions                           {GAUSS}
%  .ft      Tube length [m]                           {164mm}
%  .Mt      Transversal magnification                 {60}
%  .NA      Numerical aperture                        {1.40}
%  .Na      Samples across aperture                   {0:auto}
%  .nl      Refraction index of the lens              {1.52}
%  .nm      Refraction index of the immersion         {1.52}
%  .ns      Refraction index of the sample            {1.52}
%  .lo      Wavelength in free space [m]              {642nm}
%  .rx                                                {3um}
%  .ry      Region half axes [m]                      {3um}
%  .rz                                                {5um}
%  .xo
%  .yo      Region center [m]                         {origin}
%  .zo
%
% out    Result data
%  .dr      Radial sampling [m]                       {20nm}
%  .dz      Axial sampling [m]                        {50nm}
%
%  .Po      Laser power [W]                           {1mW}
%  .wa      Beam waist at aperture [m]                {5mm}
%
%The field functions define the amplitude and phase of the
%incident, respectively transmitted field. The following
%functions are already defined. Users can add customized
%modifiers if respecting the default arguments.
%
%     function E=custom(sys,r,t,p)
%
%      sys    System data (see examples below)
%      E      Electric field [V/m]
%      r      Radial position [Da/2]
%      t      Incidence angle [rad]
%      p      Polar angle [rad]
%
%Each field function can modify the amplitude, phase and
%polarization of the incident field. For best performance,
%polarization modifiers should be called last. The default
%field has unit amplitude, zero phase and linear x-polari-
%zation.
%
%
%CIRCULAR: Circular polarization.
%
%EXTERN:  Import the field from a MATLAB file or function.
%         Passing "Extern.mat" as file imports the default
%         MATLAB image.
%
%  .file    MATLAB file containing E, x and y
%  .func    MATLAB function returning E(x,y)
%
%        The coordinates are given with respect
%        to the aperture radius Ra=Da/2.
%
%GAUSS:   2D-Gaussian amplitude.
%
%  .Po      Laser power [W]                           {1mW}
%  .wa      Beam waist at aperture [m]                {5mm}
%
%LAMBDA4: Lambda/4 plate.
%
%  .p4      Fast axis angle [rad]
%
%LINEAR:  Linear polarization.
%
%  .pl      Polarizer angle [rad]
%
%OFFSET:  Lateral offset of next element (misalignment).
%
%  .ox      Lateral offset(s) [m]
%  .oy
%
%POLYGON: Polygon aperture.
%
%  .Nc      Number of polygon corners
%  .Rc      Radius at polygon corners [m]
%
%SEIDEL:  Aberrated wavefront (Seidel coefficients).
%
%  .Wr      Wavefront aberrations [wavelengths at aperture edge]
%    .astig    astigmatism
%    .coma     coma
%    .curv     field curvature
%    .defoc    defocus
%    .dist     distortion
%    .pist     piston
%    .spher    spherical
%    .tiltx    tilt x
%    .tilty    tilt y
%
%STRAT:   Transmission through stratified media. The wavefront
%         is returned with respect to the media the objective
%         was designed for. The media are inserted between the
%         immersion and the sample. The design focus is at z=0
%         (z pointing away from the objective).
%
%  .eff     Effective situation
%    .nl       Refraction indices of the layers
%    .hl       Thickness of the layers [m]
%  .ref     Reference/design situation
%    .nl       Refraction indices of the layers
%    .hl       Thickness of the layers [m]

%Copyright ? Marcel Leutenegger, 2003-2007, ?cole Polytechnique F?d?rale de Lausanne (EPFL),
%Laboratoire d'Optique Biom?dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
%
%    This library is free software; you can redistribute it and/or modify it under
%    the terms of the GNU Lesser General Public License as published by the Free
%    Software Foundation; version 2.1 of the License.
%
%    This library is distributed in the hope that it will be useful, but WITHOUT ANY
%    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License along
%    with this library; if not, write to the Free Software Foundation, Inc.,
%    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
%
function [sys,out,opt]=effInit(s,o,opt)
if nargin < 1 | ~isstruct(s)
   s=struct([]);
end

sys.Da=option(s,'Da',6.5e-3);
sys=optional(sys,s,'Dc');
sys=optional(sys,s,'eff');
sys.Ei=option(s,'Ei',{'gauss', 'circular'});
sys=optional(sys,s,'file');
sys.ft=option(s,'ft',250e-3);   
sys=optional(sys,s,'func');
sys.lo=option(s,'lo',635e-9); % Wavelength for excitation. 
sys.loem=option(s,'lo',700e-9); % Wavelength for emission.
sys.Mt=option(s,'Mt',60);
sys.NA=option(s,'NA',1.4);
sys=optional(sys,s,'Na');
sys=optional(sys,s,'Nc');
sys.nl=option(s,'nl',1.518);
sys.nm=option(s,'nm',1.518);
sys.ns=option(s,'ns',1.518);
sys=optional(sys,s,'ns');
sys=optional(sys,s,'p4');
sys.pl = opt.polAngle;
sys.Po=option(s,'Po',1e-5); 
sys=optional(sys,s,'Rc');
sys=optional(sys,s,'ref');
sys.rx=option(s,'rx',3e-6);
sys.ry=option(s,'ry',3e-6);
sys.rz=option(s,'rz',5e-6);
sys.rx = opt.radiusCanvas;  
sys.ry = opt.radiusCanvas;
sys.rz = opt.depthCanvas;
% sys.rx = radiusCanvas;  
% sys.ry = radiusCanvas;
% sys.rz = depthCanvas;

sys=optional(sys,s,'xo');
sys=optional(sys,s,'yo');
sys=optional(sys,s,'zo');
sys.wa=option(s,'wa',7e-3);

sys=optional(sys,s,'Wr');
% Parameters for Stratified layers. 
    sys.eff.nl = 1.4;  % Refractive index of stratified layer
    sys.eff.hl = 1e-6; % Thickness of the inserted layer [m]
    sys.ref.nl = 1.518;  % Refractive index of the reference layer
    sys.ref.hl = 0e-6; % Thickness of the inserted layer [m]
% Tilt parameters for Seidel
sys.Wr.tiltx = option(sys,'tiltx',0);
sys.Wr.tilty = option(sys,'tilty',0);
sys.W.coma = option(sys,'coma',0);
sys.Wr.defoc=option(sys,'defoc',0);
sys.Wr.spher=option(sys,'spher',0); 
sys.Wr.astig=option(sys,'astig',0);
sys.Wr.dist=option(sys,'dist',0);

%%
sys.xo=0;
sys.yo=0;
sys.zo=0;
sys.ox= 0;
sys.oy= 0;



%%
if ischar(sys.Ei)
   sys.Ei={sys.Ei};
end
if nargout < 2
   return;
end
%
% Simulation data
%
if nargin < 2 | ~isstruct(o)
   o=struct([]);
end
% out.dr = pixSize;
% out.dz = pixSize;
out.dr = opt.pixSize;
out.dz = opt.pixSize;


