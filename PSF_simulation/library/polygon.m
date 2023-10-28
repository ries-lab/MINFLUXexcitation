%Polygon aperture.
%
% sys    System data
%  .Da      Aperture diameter [m]
%  .Nc      Number of polygon corners
%  .Rc      Radius at polygon corners [m]
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
function E=polygon(sys,E,r,t,p)
t=2*pi/sys.Nc;
R=2*sys.Rc/sys.Da;
if sys.Nc > 2
   R=R*sqrt(0.5 + 0.5*cos(t));
end
e=pow2(1 + tanh(2*length(r)*(R - r.*cos(p))),-sys.Nc);
for n=1:sys.Nc-1
   e=e.*(1 + tanh(2*length(r)*(R - r.*cos(p + n*t))));
end
E=E.*repmat(e,size(E,1),1);
