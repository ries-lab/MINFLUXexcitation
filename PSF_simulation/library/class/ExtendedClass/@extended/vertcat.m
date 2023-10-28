%Vertical concatenation.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=vertcat(varargin)
o=cat(1,varargin{:});
