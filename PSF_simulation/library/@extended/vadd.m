%o=vadd(s,t)
%-----------
%
%Vector addition.
%
%  Leutenegger Marcel © 19.7.2005
%
%Input:
% s      vector
% t      vector
%
%Output:
% o      vector: s+t
%
function o=vadd(s,t)
switch nargin
case 0
   fprintf('\nVector addition.\n\n\tLeutenegger Marcel © 19.7.2005\n');
case 2
   if isempty(s) | isempty(t)
      o=[];
   else
      d=size(s);
      e=size(t);
      if d(1) ~= 3 | e(1) ~= 3
         error('Incompatible dimensions.');
      end
      if isequal(d,e)
         o=s+t;
      elseif numel(s) == 3
         o=reshape([s(1)+t(1,:);s(2)+t(2,:);s(3)+t(3,:)],e);
      elseif numel(t) == 3
         o=reshape([s(1,:)+t(1);s(2,:)+t(2);s(3,:)+t(3)],d);
      else
         error('Incompatible dimensions.');
      end
   end
otherwise
   error('Incorrect number of arguments.');
end
