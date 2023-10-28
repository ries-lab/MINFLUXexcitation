%effPlotQuadrants(out,cut,x,y,z)
%-------------------------------
%
%Plot the intensity near the focus of an objective in 7 of 8 quadrants.
%
% out    Result data
%  .dr      Sampling step [m]
%  .x
%  .y       Coordinates [m]
%  .z
%  .E       Electric field [V/m]
%
% cut    Display ratio(s) I/max(I)  {exp(-1:-1:-4)}
% x
% y      Center coordinates [m]     {0}
% z

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
function effPlotQuadrants(out,cut,x,y,z)
if nargin < 2 | isempty(cut)    %Determine whether array is empty.
   cut=exp(-1:-1:-4);
end
if nargin < 5   %Number of arguments.
   x=0;
   y=0;
   z=0;
end
if ~isfield(out,'I')    %Determine whether input is structure array field.
   out.I=vnorm(out.E);
end
out.E=[];
out.Et=[];                       % reduce memory load
out=reduce(out,min(cut));
out.I=permute(out.I,[3 2 4 1]);
out=precision(out,0);            % convert to double
out.x=1e6*out.x;
out.y=1e6*out.y;                 % coordinates in um
out.z=1e6*out.z;
x=center(out.x,1e6*x);
y=center(out.y,1e6*y);           % center sample
z=center(out.z,1e6*z);
Ix=squeeze(out.I(out.y >= y,out.x == x,out.z >= z));
Iy=squeeze(out.I(out.y == y,out.x >= x,out.z >= z));
Iz=squeeze(out.I(out.y >= y,out.x >= x,out.z == z));
[X,Y,Z]=meshgrid(out.x,out.y,out.z);
I=out.I;
I(out.y > y,out.x > x,out.z > z)=nan;
x=out.x(out.x >= x);
y=out.y(out.y >= y);             % coordinates on sections
z=out.z(out.z >= z);
f=figure('NumberTitle','off','Name','Intensity');
a=axes('Parent',f,'DataAspectRatio',[1 1 1]);
view(130,25);
for k=1:numel(cut)  %This returns the number of elements, n, in array A.
   c=[1 (k-1)/max(1,numel(cut)-0.99) 0];
   isonormals(X,Y,Z,out.I,patch(isosurface(X,Y,Z,I,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',c,'Parent',a));   %Compute normals of isosurface vertices(Cho-ten or Ko-ten).
   p=contourc(x,y,Iz,cut([k k]));
   while ~isempty(p)
      m=p(2,1);
      n=2:m+1;
      line(p(1,n),p(2,n),repmat(z(1),1,m),'Color',c,'Parent',a);
      p=p(:,2+m:end);
   end
   %Below script was originally active.
   p=contourc(z,x,Iy,cut([k k]));
   
   while ~isempty(p)    %Determine whether array is empty.
      m=p(2,1);
      n=2:m+1;
      line(p(2,n),repmat(y(1),1,m),p(1,n),'Color',c,'Parent',a);
      p=p(:,2+m:end);
   end
   %Below script was originally active.
   p=contourc(z,y,Ix,cut([k k]));
   while ~isempty(p)
      m=p(2,1);
      n=2:m+1;
      line(repmat(x(1),1,m),p(2,n),p(1,n),'Color',c,'Parent',a);
      p=p(:,2+m:end);
   end
end
axis('tight');
lighting('phong');
grid('on');
camlight;


%Determine the center sample position.
%
% r      Positions [m]
% c      Center sample [m]
%
function c=center(r,c)
n=sum(r <= c);
if n
   c=r(n);
end
