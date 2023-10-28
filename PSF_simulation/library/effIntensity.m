%out=effIntensity(sys,out)
%-------------------------
%
%Electromagnetic intensity in the focus region.
%
%Input:
% sys    System data
%  .nm      Refraction index of the immersion
%  .ns      Refraction index of the sample
% out    Result data
%  .E       Electric field [V/m]
%
%Output:
% out    Result data
%  .I       Intensity [W/m^2]

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
function out=effIntensity(sys,out)
if isfield(sys,'ns') & sys.ns < sys.nm
   n=sys.ns;
else
   n=sys.nm;
end
out.I=sqrt(e0/u0)*real(n)/2*vnorm(out.E);
