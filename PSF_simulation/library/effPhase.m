%out=effPhase(out,sys)
%---------------------
%
%Phase map of the focus field. Unwraps the phase of the
%focus field by trying to minimize phase jumps greater
%than pi between adjacent points.
%
%Input:
% sys    System data
%  .xo
%  .yo      Region center [m]
%  .zo
% out    Result data
%  .x
%  .y       Coordinates [m]
%  .z
%  .E       Electric field [V/m]
%
%Output:
% out    Result data
%  .p       Phase map [rad]

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
function out=effPhase(out,sys)
m=find(out.x == option(sys,'xo',0));    %Gives an error "Reference to non-existent field 'x'.".
n=find(out.y == option(sys,'yo',0));   % find reference point
k=find(out.z == option(sys,'zo',0));
[p,M,N,K]=size(out.E);                 % number of points
p=angle(out.E);

z=unwrap(p(:,m,n,:),[],4);             % reference axis
z=z - repmat(z(:,1,1,k),[1 1 1 K]);
x=unwrap(p(:,:,n,:),[],2);
x=x - repmat(x(:,m,1,:) - z,[1 M 1 1]);

p=unwrap(p,[],3);
p=p - repmat(p(:,:,n,:) - x,[1 1 N 1]);
if K > 1
   x=unwrap(p(:,:,:,[k+1 k]),[],4);    % correct focus plane
   p(:,:,:,k)=x(:,:,:,2);
end
clear k m n x z K M N;
out.p=p;
