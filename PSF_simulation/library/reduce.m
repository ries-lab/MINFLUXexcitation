%Reduce region to necessary range for plotting.
%
% out    Simulation results
% cut    Minimum intensity

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
function out=reduce(out,cut)
I=max(out.I(:));
z=uint8(out.I >= I*cut);
x=points(z,3,4);
y=points(z,2,4);
z=points(z,2,3);
out.x=out.x(x);
out.y=out.y(y);
out.z=out.z(z);
out.I=out.I(:,x,y,z)/I;


%Test which grid points are necessary.
%
% o      Logical array representing necessary points
% m/n    Dimensions to reduce (m/n = {1:x,2:y,3:z})
%
function o=points(o,m,n)
o=any(any(o,m),n);
o=o(:);
m=o([2:end 1]);
m(end)=0;
n=o([1 1:end-1]);
n(1)=0;
o=o | m | n;
