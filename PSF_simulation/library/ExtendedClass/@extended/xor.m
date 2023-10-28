%Logical exclusive or.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=xor(s,t)
o=(s ~= 0) ~= (t ~= 0);
