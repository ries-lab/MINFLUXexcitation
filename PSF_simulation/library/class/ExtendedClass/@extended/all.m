%True if all elements of a vector are nonzero.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=all(s,d)
if nargin < 2
   o=all(s ~= 0);
else
   o=all(s ~= 0,double(d));
end
