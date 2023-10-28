%Number of dimensions.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=ndims(s)
o=max(2,ndims(s.extended)-1);
