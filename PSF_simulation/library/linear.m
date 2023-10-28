%Linear polarization.
%
% sys    System data
%  .pl      Polarizer angle [rad]
%
% E      Electric field [V/m]
% r      Radial position [Da/2]
% t      Incidence angle [rad]
% p      Polar angle [rad]

%Copyright � Marcel Leutenegger, 2003-2007, �cole Polytechnique F�d�rale de Lausanne (EPFL),
%Laboratoire d'Optique Biom�dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
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
function E=linear(sys,E,r,t,p)
c=cos(sys.pl);
s=sin(sys.pl);
if size(E,1) > 1
   E=[c.*c.*E(1,:,:) + c.*s.*E(2,:,:);
      c.*s.*E(1,:,:) + s.*s.*E(2,:,:);
      E(3,:,:)];
else
   E=[c.*E;s.*E;zeros(size(E))];
end
end
