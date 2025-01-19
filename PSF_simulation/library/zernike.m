% Jonas Ries: zernike based on Seidel function
% 
%Aberrated wavefront (Seidel coefficients).
%
% sys    System data
%  .Da      Aperture diameter [m]
%  .Zr      Zernike polynomials; [n,l,amplitude;...] n: radial (0...) l:
%  azimuthal. In units of Pi. Or better in rad?
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
function E=zernike(sys,E,r,t,p)
% e=0;
% t=r.*r;
Zr=sys.Zr;

r(r>1)=1; %calculated only for r<1, then extended 
e=zeros(size(r));
for k=1:size(Zr,1)
    zpol=zernfun(Zr(k,1),Zr(k,2),r,p);
    zpolr=reshape(zpol,size(r));
    e=e+zpolr*Zr(k,3);
end
%
% On axis aberrations
% %
% if isfield(W,'defoc')
%    e=W.defoc*t;                  % r^2
% end

E=E.*repmat(exp(2i*pi*e),size(E,1),1);  %B = repmat(A,m,n) returns an array containing m x n copies of A in the row and column dimensions.
