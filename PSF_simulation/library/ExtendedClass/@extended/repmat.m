%Replicate and tile an array.
%
%    Marcel Leutenegger © 23.5.2005
%
function o=repmat(o,m,n)
m=double(m);
if nargin < 3
   if numel(m) > 1
      n=1;
   else
      n=m;
   end
else
   n=double(n);
end
if isa(o,'extended')
   o.extended=repmat(o.extended,[1 m n]);
else
   o=repmat(o,[m n]);
end
