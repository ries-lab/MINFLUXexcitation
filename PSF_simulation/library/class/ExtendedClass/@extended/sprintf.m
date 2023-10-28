%Write formatted data to string.
%
%	Marcel Leutenegger © 2.5.2005
%
function [s,m]=sprintf(varargin)
for n=2:nargin
   if isa(varargin{n},'extended')
      varargin{n}=double(varargin{n});
   end
end
if nargout < 2
   s=sprintf(varargin{:});
else
   [s,m]=sprintf(varargin{:});
end
