%Logical or.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=or(s,t)
o=(s ~= 0) | (t ~= 0);
