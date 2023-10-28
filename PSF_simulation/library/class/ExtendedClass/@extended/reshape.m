%Change size.
%
%	Marcel Leutenegger © 26.12.2007
%
function o=reshape(o,varargin)
for n=1:nargin-1
   varargin{n}=double(varargin{n});
end
if isa(o,'extended')
   if nargin > 2
      o.extended=reshape(o.extended,10,varargin{:});
   else
      o.extended=reshape(o.extended,[10 varargin{1}]);
   end
else
   o=reshape(o,varargin{:});
end
