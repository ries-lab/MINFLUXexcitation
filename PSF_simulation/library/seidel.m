%Aberrated wavefront (Seidel coefficients).
%
% sys    System data
%  .Da      Aperture diameter [m]
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
%  .nm      Refraction index of the medium
%  .lo      Wavelength in free space [m]
%
% E      Electric field [V/m]
% r      Radial position [Da/2]
% t      Incidence angle [rad]
% p      Polar angle [rad]

%Copyright © Marcel Leutenegger, 2003-2007, École Polytechnique Fédérale de Lausanne (EPFL),
%Laboratoire d'Optique Biomédicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
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
function E=seidel(sys,E,r,t,p)
e=0;
t=r.*r;
W=sys.Wr;
%
% On axis aberrations
%
if isfield(W,'defoc')
   e=W.defoc*t;                  % r^2
end
if isfield(W,'spher')
   e=e + W.spher*t.*t;           % r^4
end
%
% Off axis aberrations (q: tilt)
%
x=option(W,'tiltx',0);
y=option(W,'tilty',0);
q=x.*x + y.*y;
if q 
   r=r.*cos(p - atan2(y,x));    %returns the four-quadrant inverse tangent (tan-1) of Y and X.atan2(4,-3) = 2.2143

   p=sqrt(q);
   e=e + p.*r;                   % q*r*cos(p)   Tip/Tilt
   if isfield(W,'astig')
      e=e + W.astig*q.*r.*r;     % q^2*r^2*cos(p)^2
   end
   if isfield(W,'coma')
      e=e + W.coma*p.*r.*t;      % q*r^3*cos(p)
   end
   if isfield(W,'curv')
      e=e + W.curv*q.*t;         % q^2*r^2
   end
   if isfield(W,'dist')
      e=e + W.dist*p.*q.*r;      % q^3*r*cos(p)
   end
   if isfield(W,'pist')
      e=e + W.pist*q;            % q^2
   end
end
E=E.*repmat(exp(2i*pi*e),size(E,1),1);  %B = repmat(A,m,n) returns an array containing m x n copies of A in the row and column dimensions.
