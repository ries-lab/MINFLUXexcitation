%effPlotIntensity(out,cut)
%-------------------------
%
%Plot the intensity near the focus of an objective.
%
% out    Result data
%  .dr      Sampling step [m]
%  .x
%  .y       Coordinates [m]
%  .z
%  .E       Electric field [V/m]
%
% cut    Display ratio(s) I/max(I) {exp(-1:-1:-4)}

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
function effPlotIntensity(out,cut)
if nargin < 2 | isempty(cut)
   cut=exp(-1:-1:-4);
end
if ~isfield(out,'I')
   out.I=vnorm(out.E);
end
out.E=0;
out=reduce(out,min(cut));
out.I=permute(out.I,[3 2 4 1]);
out=precision(out,0);
[X,Y,Z]=meshgrid(out.x*1e6,out.y*1e6,out.z*1e6);
f=figure('NumberTitle','off','Name','Intensity');
a=axes('Parent',f,'DataAspectRatio',[1 1 1]);
view(35,45);
for k=1:numel(cut)
   isonormals(X,Y,Z,out.I,patch(isosurface(X,Y,Z,out.I,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
end
axis('tight');
lighting('phong');
grid('on');
camlight;
