%o=precision(o,p)
%----------------
%
%Set precision of all numerical parameters in a structure.
%
%Input:
% o      Parameter structure
% p      Precision class or constructor
%
%Output:
% o      Converted parameter structure
%

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
function o=precision(o,p)
if nargin > 1
   p=['=' class(p) '(v);'];
   if isstruct(o)
      n=fieldnames(o);
      for m=1:numel(n)
         v=getfield(o,n{m});
         if isnumeric(v)
            eval(['o.' n{m} p]);
         end
      end
   else
      v=o;
      eval(['o' p]);
   end
end
