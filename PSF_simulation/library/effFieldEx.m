%[out,err]=effFieldEx(sys,out,opt)
%---------------------------------
%
%Electric field near the focus of an objective with the
%chirp z-transform implementation of the extended Debye
%integral.
%
%Input:
% sys    System data
% out    Result data
% opt    Options
%  .Et      Display transmitted field           {0}
%  .Ef      Display focus field                 {1}
%  .mem     System memory (RAM) [MB]        {512MB}
%
%Output:
% out    Result data
%  .x
%  .y       Coordinates [m]
%  .z
%  .E       Electric field [V/m]
% err    Error message

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
function [out,err]=effFieldEx(sys,out,m)
if nargin < 3 | ~isstruct(m)
   m=struct([]);
end
opt.Et=option(m,'Et',0);
opt.Ef=option(m,'Ef',1);
opt.mem=option(m,'mem',512);
m=ceil(sys.rz/out.dz);
if isfield(sys,'zo')
   z=sys.zo + out.dz*(-m:m);
else
   z=out.dz*(0:m);
   sys.zo=0;
end
if ~isfield(sys,'ns')
   sys.ns=sys.nm;
elseif sys.ns < sys.nm
   z=z(z >= 0);
   if isempty(z)
      error('Nothing to do: region is only in immersion.');
   end
end
nz=numel(z);
z=reshape(z,[1 1 nz]);
[Et,kz,m]=field(sys,out);
k=pi*sys.ns/sys.lo;
dko=2*pi/m;
out.Et=Et;
if opt.Et
   imagesp(vabs(Et),'Field Et');
end
m=ceil(sys.rx/out.dr);
if isfield(sys,'xo')
   out.x=sys.xo + out.dr*(-m:m).';
else
   out.x=out.dr*(0:m).';
end
m=ceil(sys.ry/out.dr);
if isfield(sys,'yo')
   out.y=sys.yo + out.dr*(-m:m);
else
   out.y=out.dr*(0:m);
end
nx=numel(out.x);
ny=numel(out.y);
out.z=z;
if isa(Et,'single')
   s=24/1024^2;
else
   s=48/1024^2;
end
if s*nx*ny*nz > max(opt.mem-128,opt.mem/2)
   error('Out of memory.');
end
%
% Prepare chirp z-transform: constants
%
[s,m,n]=size(Et);
M=optdim(m+nx-1);
N=optdim(n+ny-1);
f=sys.nm*sys.ft/sys.Mt;
P=reshape(f.^2 + repmat(out.x.^2,1,ny) + repmat(out.y.^2,nx,1),[1 nx ny]);
s=f./(f + sys.ns/sys.nm*out.z);
%
% Compute field with chirp z-transform
%
out.E=Et(1);
out.E(3,nx,ny,nz)=0;
h=statusbar('Electric field ...');
for z=1:nz
   dk=s(z)*dko;
   E=ifft(vmul(fft(vmul(Et,exp(i*s(z)*out.z(z)*kz).*repmat(cis(-dk/2*(0:m-1).^2 - dk*out.x(1)/out.dr*(0:m-1)),[1 1 n])),M,2),repmat(fft(cis(dk/2*(1-m:nx-1).^2),M),[1 1 n])),M,2);
   E=ifft(vmul(fft(vmul(E(:,m:m+nx-1,:),reshape(cis(dk/2*((m-1)*(0:nx-1) - (0:nx-1).^2).')*cis(-dk/2*(0:n-1).^2 - dk*out.y(1)/out.dr*(0:n-1)),[1 nx n])),N,3),repmat(reshape(fft(cis(dk/2*(1-n:ny-1).^2),N),[1 1 N]),1,nx)),N,3);
   E=vmul(E(:,:,n:n+ny-1),s(z)*cis(k*(out.z(z) + (P + out.z(z).^2)/(f + out.z(z)))).*repmat(reshape(cis(dk/2*((n-1)*(0:ny-1) - (0:ny-1).^2)),[1 1 ny]),1,nx));
   out.E(:,:,:,z)=E;
   if opt.Ef
      imagesp(vabs(E),'Field Ef');
   end
   clear('E');
   if isempty(statusbar(z/nz,h))
      err='User abort.';
      return;
   end
end
delete(h);
err='';
