%effPlotSections(out,x,y,z)
%--------------------------
%
%Plot cross-sections (principal planes) through the focus.
%
% out    Result data
%  .dr      Radial sampling step [m]
%  .dz      Axial sampling step [m]
%  .x
%  .y       Coordinates [m]
%  .z
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
function effPlotSections(out,x,y,z)
out.I=squeeze(out.I);
if nargin < 2
   [x,y,z]=ind2sub(size(out.I),find(out.I == max(out.I(:))));
   [a,f]=min(abs(out.x(x))); x=x(f);
   [a,f]=min(abs(out.y(y))); y=y(f);
   [a,f]=min(abs(out.z(z))); z=z(f);
else
   [a,x]=min(abs(out.x-x));
   [a,y]=min(abs(out.y-y));
   [a,z]=min(abs(out.z-z));
end
Ix=log(squeeze(out.I(x,:,:)));
Iy=log(squeeze(out.I(:,y,:)));
Iz=log(squeeze(out.I(:,:,z)));
m=min([Ix(:);Iy(:);Iz(:)]);
n=max([Ix(:);Iy(:);Iz(:)]);
x=out.x*1e6;
y=out.y*1e6;
z=out.z*1e6;

[X,Y,Z]=size(out.I);
K=[x(1)-out.dr/2 x(end)+out.dr/2];
M=[y(1)-out.dr/2 y(end)+out.dr/2];
N=[z(1)-out.dz/2 z(end)+out.dz/2];

a.CLim=[m n];       %Gives an error "Field assignment to a non-structure array object".
a.FontName='Times New Roman';
a.Parent=figure('NumberTitle','off','Name','Cross-sections','Color','white','Colormap',jetsp);
a.TickDir='out';
a.Units='pixel';
a.YDir='normal';

h=axes(a,'Position',[50 50 X Y]);  % Original
%h=axes('Position',[50 50 X Y]);
%imagesc(x,y,Iz.','Parent',h);
imagesc(x,y,Iz.');
set(h,a);
xlabel('x');
ylabel('y');

%h=axes(a,'Position',[70+X 50 Z Y]);
h=axes('Position',[70+X 50 Z Y]);
%imagesc(z,y,Ix,'Parent',h);
imagesc(z,y,Ix);
set(h,a,'YTickLabel',{''});
xlabel('z');

%h=axes(a,'Position',[50 70+Y X Z]);
h=axes('Position',[50 70+Y X Z]);
%imagesc(x,z,Iy.','Parent',h);
imagesc(x,z,Iy.');
set(h,a,'XTickLabel',{''});
ylabel('z');

h=axes(a,'Position',[100+X+Z 50 20 20+Y+Z]);
surf([0 1],double(exp([m n])),zeros(2),double([m m;n n]),'Parent',h,'FaceColor','interp','EdgeColor','black');
set(h,a,'Box','on','XTick',[],'YAxisLocation','right','YScale','log','YLim',double(exp([m n])));
ylabel('I');
view(0,90);
