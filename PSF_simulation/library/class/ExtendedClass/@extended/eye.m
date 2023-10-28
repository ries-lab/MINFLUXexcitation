%Identity matrix.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=eye(varargin)
o=zeros(varargin{:});
t=true(o);
t=tril(t) & triu(t);
v=ones(extended(min(size(t))),1);
o.extended(:,t)=v.extended;
