%True if arrays are numerically equal.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=isequal(s,varargin)
switch nargin
case 0
   fprintf('\nTrue if arrays are numerically equal.\n\n\tMarcel Leutenegger © 2.5.2005\n\n');
case 1
   error('Wrong number of input arguments');
otherwise
   if ~isa(s,'extended')
      s=extended(s);
   end
   d=size(s);
   for n=1:nargin-1
      o=isequal(d,size(varargin{n}));
      if o
         o=s == varargin{n};
         o=all(o(:));
      end
      if ~o
         break;
      end
   end
end
