%Write formatted data to file.
%
%	Marcel Leutenegger © 2.5.2005
%
function c=fprintf(f,varargin)
if ischar(f)
   n=1;
else
   n=2;
end
for n=n:nargin-1
   if isa(varargin{n},'extended')
      varargin{n}=double(varargin{n});
   end
end
if nargout
   c=fprintf(f,varargin{:});
else
   fprintf(f,varargin{:});
end
