%Matrix power.
%
%  Marcel Leutenegger © 25.4.2007
%
function o=mpower(s,t)
if isempty(s) | isempty(t)
   o=[];
elseif numel(s)*numel(t) == 1
   o=s.^t;
else
   o=extended(double(s)^double(t));
end
