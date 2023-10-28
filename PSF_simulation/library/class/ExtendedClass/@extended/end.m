%Index of last element in a dimension.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=end(s,k,n)
m=ndims(s.extended);
k=k+1;
n=n+1;
if k > m
   o=1;
else
   if k < n
      o=size(s.extended,k);
   else
      o=size(s.extended);
      o=prod(o(k:m));
   end
end
