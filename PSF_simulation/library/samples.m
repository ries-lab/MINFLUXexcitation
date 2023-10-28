%Sampling over aperture and lateral wavevector.
%
% sys    System data
% out    Simulation results
%
% M      Samples on aperture radius
% N      Samples on lateral wavevector

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
function [M,N]=samples(sys,out)
if sys.ns < sys.nm
   M=max(100,21*sys.NA*sys.ns*(sys.zo + sys.rz)^2/sys.lo^2);
else
   M=max(50,7*sys.NA^2/sqrt(sys.ns^2 - sys.NA^2)*(abs(sys.zo) + sys.rz)/sys.lo);
end
if isfield(sys,'da')
   M=max(M,sys.Da/sys.da);
end
if isfield(sys,'Na')
   M=max(M,sys.Na);
end
N=ceil(M*max(2,sys.lo/2/sys.NA/out.dr));
M=N*sys.NA*out.dr/sys.lo;
