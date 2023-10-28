%Concatenation.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=cat(d,varargin)
if isa(d,'extended')
   o=cat(double(d),varargin{:});
else
   data=cell(1,nargin-1);
   for n=1:nargin-1
      if isa(varargin{n},'extended');
         data{n}=varargin{n}.extended;
      else
         data{n}=extended(varargin{n});
         data{n}=data{n}.extended;
      end
   end
   o=extended;
   o.extended=cat(max(2,d+1),data{:});
end
