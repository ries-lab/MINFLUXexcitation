%Right matrix divide.
%
%    Marcel Leutenegger © 23.5.2005
%
function o=mrdivide(s,t)
if isempty(s)
   o=[];
else
   if numel(t) == 1
      o=s./t;
   else
      o=single(double(s)/double(t));
   end
end
