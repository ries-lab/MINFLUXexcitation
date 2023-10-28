%Horizontal concatenation.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=horzcat(varargin)
o=cat(2,varargin{:});
