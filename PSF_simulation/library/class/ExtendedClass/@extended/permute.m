%Permute array dimensions.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=permute(o,d)
d=double(d);
if isa(o,'extended')
   o.extended=permute(o.extended,[1 1+d]);
else
   o=permute(o,d);
end
