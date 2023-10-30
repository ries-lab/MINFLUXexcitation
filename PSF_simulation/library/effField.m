%[out,err]=effField(sys,out,opt)
%-------------------------------
%
%Electric field near the focus of an objective with the
%chirp z-transform implementation of the Debye integral.
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

%Copyright � Marcel Leutenegger, 2003-2007, �cole Polytechnique F�d�rale de Lausanne (EPFL),
%Laboratoire d'Optique Biom�dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
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
%function [out,err]=effField(sys,out,m)
function [out,err]=effField(sys,out,opt)
% if nargin < 3 | ~isstruct(m)
%    m=struct([]);
if nargin < 3 | ~isstruct(opt)
   opt=struct([]);
end
% opt.Et=option(m,'Et',0);    % Display Et?
% opt.Ef=option(m,'Ef',0);    % Display Ef?
% opt.mem=option(m,'mem',4096);

m=ceil(sys.rz/out.dz);  % ceil rounds up to integer.
if isfield(sys,'zo')
   z=sys.zo + out.dz*(-m:m);    % Setting the region center.
else
   z=out.dz*(0:m);
   sys.zo=0;
end


if ~isfield(sys,'ns')   % if not the sample RI, then...
   sys.ns=sys.nm;           % Set the RI of sample to the same of immersion.
elseif sys.ns < sys.nm
   z=z(z >= 0);
   if isempty(z)
      error('Nothing to do: region is only in immersion.');
   end
end


nz=numel(z);
z=reshape(z,[1 1 nz]);
[Et,kz,m]=field(sys,out);   % Function to calculate the transmission through objective lens.
dk=2*pi/m;
out.Et=Et;

if opt.Et
   imagesp(vabs(Et),'Field Et');    % Displaying transmission field.
end

if z(1)
   Et=vmul(Et,exp(i*z(1)*kz));
end

m=ceil(sys.rx/out.dr);
if isfield(sys,'xo')
   out.x=sys.xo + out.dr*(-m:m).';      % Setting the region center in X direction.
else
   out.x=out.dr*(0:m).';
end


m=ceil(sys.ry/out.dr);
if isfield(sys,'yo')
   out.y=sys.yo + out.dr*(-m:m);        % Setting the region center in Y direction.
else
   out.y=out.dr*(0:m);
end


kz=exp(i*out.dz*kz);
nx=numel(out.x);
ny=numel(out.y);
out.z=z;
if isa(Et,'single')
   s=24/1024^2;     % Is this a pixel size??
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
Et=vmul(Et,repmat(cis(-dk/2*(0:m-1).^2 - dk*out.x(1)/out.dr*(0:m-1)),[1 1 n]));
a=reshape(cis(dk/2*((m-1)*(0:nx-1) - (0:nx-1).^2).')*cis(-dk/2*(0:n-1).^2 - dk*out.y(1)/out.dr*(0:n-1)),[1 nx n]);
b=repmat(reshape(cis(dk/2*((n-1)*(0:ny-1) - (0:ny-1).^2)),[1 1 ny]),1,nx);
v=repmat(fft(cis(dk/2*(1-m:nx-1).^2),M),[1 1 n]);
w=repmat(reshape(fft(cis(dk/2*(1-n:ny-1).^2),N),[1 1 N]),1,nx);
m=m:m+nx-1;
n=n:n+ny-1;
%
% Compute field with chirp z-transform
%
out.E=Et(1);
out.E(3,nx,ny,nz)=0;
if opt.calbar
f=statusbar('Electric field ...');
end
for z=1:nz
   E=ifft(vmul(fft(Et,M,2),v),M,2);
   E=ifft(vmul(fft(vmul(E(:,m,:),a),N,3),w),N,3);
   E=vmul(E(:,:,n),b);
   out.E(:,:,:,z)=E;
   if opt.Ef            % Displaying electric field.
      E=vabs(E);
      imagesp(E,'Field Ef');
   end
   clear('E');
   if opt.calbar
   if isempty(statusbar(z/nz,f))
      err='User abort.';
      return;
   end
   end
   Et=vmul(Et,kz);
end
if opt.calbar
delete(f);
end
err='';
