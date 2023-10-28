%Right matrix divide.
%
%    Marcel Leutenegger © 23.5.2005
%
function o=mldivide(s,t)
if isempty(t)
   o=[];
else
   if numel(s) == 1
      o=s.\t;
   else
      o=single(double(s)\double(t));
   end
end
