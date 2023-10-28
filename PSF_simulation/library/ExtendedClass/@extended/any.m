%True if any element of a vector is nonzero.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=any(s,d)
if nargin < 2
   o=any(s ~= 0);
else
   o=any(s ~= 0,double(d));
end
