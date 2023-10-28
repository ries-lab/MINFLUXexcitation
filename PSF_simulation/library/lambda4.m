%Lambda/4 plate.
%
% sys    System data
%  .p4      Fast axis angle [rad]
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
function E=lambda4(sys,E,r,t,p)
c=cos(sys.p4);
s=sin(sys.p4);
if size(E,1) > 1
   E=[complex(c.*c,s.*s).*E(1,:,:) + complex(c.*s,-c.*s).*E(2,:,:);
      complex(c.*s,-c.*s).*E(1,:,:) + complex(s.*s,c.*c).*E(2,:,:);
      E(3,:,:)];
else
   E=[complex(c.*c,s.*s).*E;
      complex(c.*s,-c.*s).*E;
      zeros(size(E))];
end
