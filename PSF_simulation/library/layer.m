%Fresnel reflection and transmission coefficients at a layer.
%
% r1     Reflection coefficients
% r2
% t1     Transmission coefficients
% t2
% pz     Phase advance through layer
%
% r      Coefficients at the layer
% t

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
function [r,t]=layer(r1,t1,r2,t2,pz)
r2=r2.*pz.*pz;
t=1 + r1.*r2;
r=(r1 + r2)./t;
t=t1.*t2.*pz./t;
