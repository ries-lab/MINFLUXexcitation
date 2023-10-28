%Largest positive floating point number.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=realmax(s)
o=pow2(2-eps(s),16383);
