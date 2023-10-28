%z=cis(p,r)
%----------
%
%Sine and cosine.
%
%Input:
% p     Phase [rad]
% r     Magnitude {1}
%
%Output:
% z     r.*complex(cos(p),sin(p))
%
%Alternative:
% [c,s] To get the sine and cosine separately.

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
function [c,s]=cis(p,r)
switch nargin
case 0
   fprintf('\nSine and cosine.\n\n\tMarcel Leutenegger © 1.5.2005\n');
case 1
   c=exp(complex(0,p));
case 2
   c=r.*exp(complex(0,p));
otherwise
   error('Incorrect number of arguments.');
end
if nargout == 2
   s=imag(c);
   c=real(c);
end
