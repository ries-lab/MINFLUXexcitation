%Diagonal matrices and diagonals of a matrix.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=diag(o,k)
if nargin > 1
   k=double(k);
else
   k=0;
end
if isa(o,'extended')
   s=size(o);
   if sum(s > 1) < 2
      v=o;
      t=logical(diag(true(o),k));
      o=zeros(extended(abs(k)+prod(s)));
      o.extended(:,t)=v.extended;
   else
      t=true(o);
      t=tril(t,k) & triu(t,k);
      o.extended=o.extended(:,t);
   end
else
   o=diag(o,k);
end
