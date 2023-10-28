%Transmission through stratified media. The wavefront
%is returned with respect to the media the objective
%was designed for. The media are inserted between the
%immersion and the sample. The design focus is at z=0
%(z pointing away from the objective).
%
% sys    System data
%  .eff     Effective situation
%    .nl       Refraction indices of the layers
%    .hl       Thickness of the layers [m]
%  .ref     Reference/design situation
%    .nl       Refraction indices of the layers
%    .hl       Thickness of the layers [m]
%  .lo      Wavelength in free space [m]
%  .nm      Refraction index of the immersion
%
% t      Total deflection angle [rad]
% tp     Transmission coefficients
% ts
% b2     Propagation kz/ko

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
function [tp,ts,b2]=stratified(sys,t)
go=sys.nm*sin(t);
ko=2i*pi/sys.lo;
er=[sys.nm;sys.eff.nl(:)].^2;
[Tp,Ts]=transmission([sys.nm;sys.ref.nl(:)].^2,sys.ref.hl(:),go,ko);
[tp,ts]=transmission(er,sys.eff.hl(:),go,ko);
tp=tp.*conj(sign(Tp));
ts=ts.*conj(sign(Ts));
[Tp,Ts,b2]=interface(er(end),sys.ns.^2,go);
tp=tp.*Tp;
ts=ts.*Ts;


%Fresnel transmission coefficients at layers.
%
% er     Relative dielectric constants
% hz     Thickness of layers [m]
% go     Incidence kxy/ko
% ko     Wavevector in free space [rad/m]
%
% tp     Transmission coefficients
% ts
%
function [tp,ts]=transmission(er,hz,go,ko)
if numel(er) < 2
   tp=1;
   ts=1;
else
   n=imag(er) == 0;        % avoid undefined results in resonances
   er(n)=er(n) + 1e-9i;
   [tp,ts,bz,rp,rs]=interface(er(end-1),er(end),go);
   bz=exp(i*hz(end)*ko*bz);
   for n=numel(er)-1:-1:2
      [ta,tb,pz,ra,rb]=interface(er(n-1),er(n),go);
      pz=exp(i*hz(n-1)*ko*pz);
      [rp,tp]=layer(ra,ta,rp,tp,pz);
      [rs,ts]=layer(rb,tb,rs,ts,pz);
   end
   tp=tp.*bz;
   ts=ts.*bz;
end
