%Logical and.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=and(s,t)
o=(s ~= 0) & (t ~= 0);
