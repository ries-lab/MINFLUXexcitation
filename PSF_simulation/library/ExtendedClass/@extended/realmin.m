%Smallest positive floating point number.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=realmin(s)
o=pow2(ones(extended),-16382);
