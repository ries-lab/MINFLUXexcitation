%Execute field modifiers.
%
% sys    System parameters
% M      Samples on aperture radius
% N      Samples on lateral wavevector
%
% E      Electric field [V/m]
% r      Radial position [Da/2]
% t      Incidence angle [rad]
% p      Polar angle [rad]
% ct     cos(t)
% st     sin(t)

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
function [E,r,t,p,ct,st]=modifiers(sys,M,N)
x=sys.NA/sys.nm/M*(-ceil(M+2):ceil(M+2));
y=reshape(x,[1 size(x)]);
p=atan2(repmat(y,size(x)),repmat(x,size(y)));
r=repmat(x.*x,size(y)) + repmat(y.*y,size(x));
ct=sqrt(1 - r);
st=sqrt(r);
t=asin(st);
r=sys.nm/sys.NA*st;        % radial position in aperture
E=-0.5i*(sys.Da/2/M)^2*sys.Mt/(sys.ft*sys.lo);
clear('strat');            % stratified not yet specified
n=0;
while n < numel(sys.Ei)    % MATLAB bug: FOR does not check loop variable,
   n=n + 1;                %             so OFFSET and CENTER would fail.
   E=feval(sys.Ei{n},sys,E,r,t,p);
end
E=E.*repmat((1 + tanh(2*length(r)*(1 - r)))./ct,size(E,1),1);
