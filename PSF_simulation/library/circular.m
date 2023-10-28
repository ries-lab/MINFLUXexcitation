%Circular polarization.
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
function E=circular(sys,E,r,t,p)
if size(E,1) > 1
   E=[E(1,:,:)./sqrt(2) - i./sqrt(2).*E(2,:,:);
      i./sqrt(2).*E(1,:,:) + E(2,:,:)./sqrt(2);
      E(3,:,:)];
else
   E=[E./sqrt(2);i./sqrt(2).*E;zeros(size(E))];
end
E;
% figure;imagesc(squeeze(abs(E(1,:,:))));axis square
end