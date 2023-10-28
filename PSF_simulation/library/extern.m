%Import the field from a MATLAB file or function.
%
% sys    System data
%  .file    MATLAB file containing E, x and y
%  .func    MATLAB function returning E(x,y)
%
%        The coordinates are given with respect
%        to the aperture radius Ra=Da/2.
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
function E=extern(sys,e,r,t,p)
if isfield(sys,'file')
   load(sys.file,'E','x','y');
   E=squeeze(E);
   x=squeeze(x);
   y=squeeze(y);
   if ~isequal(size(x),size(y))
      [x,y]=ndgrid(x,y);
   end
   [p,t]=cis(squeeze(p),squeeze(r));
   E=reshape(interpn(x,y,E,p,t),size(r));
elseif isfield(sys,'func')
   [p,t]=cis(p,r);
   E=feval(sys.func,p,t);
else
   error('No external source for electric field.');
end
E(isnan(E))=0;
E=e.*repmat(E,size(e,1),1);
