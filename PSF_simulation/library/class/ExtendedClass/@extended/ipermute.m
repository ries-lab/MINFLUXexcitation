%Inverse permute array dimensions.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=ipermute(o,d)
d=double(d);
if isa(o,'extended')
   o.extended=ipermute(o.extended,[1 1+d]);
else
   o=ipermute(o,d);
end
