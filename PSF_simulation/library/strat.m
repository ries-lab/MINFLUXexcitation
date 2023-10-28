%Transmission through stratified media. The wavefront
%is returned with respect to the media the objective
%was designed for. The media are inserted between the
%immersion and the sample. The design focus is at z=0
%(z pointing away from the objective).
%
% sys    System data
%  .eff     Effective situation
%    .nl       Refraction indices of the layers
%    .hl       Thickness of the layers [m]
%  .ref     Reference/design situation
%    .nl       Refraction indices of the layers
%    .hl       Thickness of the layers [m]
%  .lo      Wavelength in free space [m]
%  .nm      Refraction index of the immersion
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
function E=strat(sys,E,r,t,p)
persistent n;
if nargin
   n(end+1)=evalin('caller','n');
else
   E=n;
end
