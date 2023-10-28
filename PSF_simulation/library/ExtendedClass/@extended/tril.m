%Lower triangular part of a matrix.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=tril(o,k)
if nargin > 1
   k=double(k);
else
   k=0;
end
o.extended(:,logical(triu(true(o),k+1)))=0;
