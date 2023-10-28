%Convert number to string.
%
%    Marcel Leutenegger © 24.5.2005
%
function t=num2str(varargin)
for n=1:nargin
   if isa(varargin{n},'extended')
      varargin{n}=double(varargin{n});
   end
end
t=num2str(varargin{:});
