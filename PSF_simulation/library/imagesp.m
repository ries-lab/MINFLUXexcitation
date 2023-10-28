%h=imagesp(i,t,c,f,a)
%--------------------
%
%Display an image with automatic (true-)color scaling.
%
% i      Image
% t      Title {Image}
% c      Corner {0:centered on screen (default), 1:upper-left,
%                2:upper-right, 3:lower-right, 4:lower-left}
% f      Figure properties
% a      Axes properties
%
% h      Axes handle

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
function h=imagesp(i,t,c,f,a)
if nargin < 2 | numel(t) < 2
   t='Image';
end
if nargin < 4 | ~isstruct(f)
   f=struct('DoubleBuffer','on','Menubar','none','NumberTitle','off','Colormap',graysp);
end
if nargin < 5 | ~isstruct(a)
   a=struct('Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren');
end
%
% Axes limits
%
i=squeeze(i);
[m,n,k]=size(i);
a.XLim=0.5+[0 m];
a.YLim=0.5+[0 n];
a.DataAspectRatio=[1 1 1];
%
%Activate or create a figure.
%
h=findobj(get(0,'Children'),'flat','Name',t);
s=get(0,'ScreenSize');
s=s(3:4)-[m n+18];
if numel(h)
   if nargin > 2 & numel(c) == 1
      switch c
      case 1
         f.Position=[4 s(2)-4 m n];
      case 2
         f.Position=[s-4 m n];
      case 3
         f.Position=[s(1)-4 33 m n];
      case 4
         f.Position=[4 33 m n];
      otherwise
         f.Position=[s/2 m n];
      end
   end
   set(h,f);
   figure(h);
   h=get(h,'Children');
   set(h,a);
else
   if nargin < 3 | numel(c) ~= 1
      c=0;
   end
   switch c
   case 1
      f.Position=[4 s(2)-4 m n];
   case 2
      f.Position=[s-4 m n];
   case 3
      f.Position=[s(1)-4 4 m n];
   case 4
      f.Position=[4 4 m n];
   otherwise
      f.Position=[s/2 m n];
   end
   h=axes('Parent',figure(f,'Name',t),a);
end
if ~isreal(i)
   i=abs(i);
end
if k == 3 & ~isa(i,'uint8')
   for k=1:3
      j=i(:,:,k);
      m=min(j(isfinite(j)));
      n=max(j(isfinite(j)));
      i(:,:,k)=255/(n-m)*(j-m);
   end
   i=uint8(i);
%    k=min(i(isfinite(i)));
%    i=uint8(255/(max(i(isfinite(i)))-k)*(i-k));
end
imagesc(permute(i,[2 1 3]),'Parent',h);
drawnow;
