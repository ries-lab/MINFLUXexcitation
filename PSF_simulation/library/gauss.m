%2D-Gaussian amplitude.
%
% sys    System data
%  .Da      Aperture diameter [m]
%  .Po      Laser power [W]
%  .wa      Beam waist at aperture [m]
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
function E=gauss(sys,E,r,t,p)
e=2*sqrt(sqrt(u0/e0)*sys.Po/pi)/sys.wa;
%e=e*exp(-((sys.Da-(0.01))/sys.wa/2).^2*r.^2);    % Changing Gaussian beam
%diameter
e=e*exp(-(sys.Da/sys.wa/2).^2*r.^2);
E=E.*repmat(e,size(E,1),1);
E;
