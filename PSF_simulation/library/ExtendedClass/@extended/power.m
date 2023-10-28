%Array power.
%
%    Marcel Leutenegger © 5.2.2005
%
function o=power(s,t)
if ~isa(s,'extended')
   s=extended(s);
end
if ~isa(t,'extended')
   t=extended(t);
end
o=exp(log(s).*t);
