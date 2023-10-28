%Electric field transmitted through an objective.
%
% sys    System data
% out    Simulation results
%
% E      Electric field [V/m]
% k      Wavevector [rad/m]
% N      Number of samples

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
function [E,k,N]=field(sys,out)
[M,N]=samples(sys,out);     %Sampling over aperture and lateral wavevector.
[E,r,t,p,ct,st]=modifiers(sys,M,N);     %Execute field modifiers.
if sys.phaseImage == true
    imtool(squeeze(angle(E(1,:,:))), [min(min(squeeze(angle(E(1,:,:))))) max(max(squeeze(angle(E(1,:,:)))))])
    imtool(squeeze(angle(E(2,:,:))), [min(min(squeeze(angle(E(2,:,:))))) max(max(squeeze(angle(E(2,:,:)))))])
end
if sys.intImage == true
    imtool(squeeze(abs(E(1,:,:))), [min(min(squeeze(abs(E(1,:,:))))) max(max(squeeze(abs(E(1,:,:)))))])
    figure;plot(abs(E(1,:,29)))
%     imtool(squeeze(real(E(1,:,:))), [min(min(squeeze(real(E(1,:,:))))) max(max(squeeze(real(E(1,:,:)))))])
    %     imtool(squeeze(abs(E(2,:,:))), [min(min(squeeze(abs(E(2,:,:))))) max(max(squeeze(abs(E(2,:,:)))))])
end
if any(strat)
   [cp,sp,k]=stratified(sys,t);         %Transmission through stratified media.
else
   [cp,sp,k]=interface(sys.nm^2,sys.ns^2,sys.NA*r);     %Fresnel reflection and transmission coefficients at an interface for the electric field (p- and s-polarisation).
end
[tp,ts]=transmission(sys,t);        %Transmission through an air-lens-air-immersion objective.
tp=tp.*cp;
ts=ts.*sp;
[cp,sp]=cis(p);     %Sine and cosine.
if size(E,1) > 1
   tp=tp.*(cp.*E(1,:,:) + sp.*E(2,:,:));
   ts=ts.*(sp.*E(1,:,:) - cp.*E(2,:,:));
else
   tp=tp.*cp.*E;
   ts=ts.*sp.*E;
end
ct=tp.*k/sys.ns;
E=cat(1,ct.*cp + sp.*ts,ct.*sp - cp.*ts,tp.*st/sys.ns);
k=2*pi/sys.lo*k;
E(isnan(E))=0;
