%Subscripted assignment.
%
%    Marcel Leutenegger © 23.5.2005
%
function o=subsasgn(o,s,v)
if ~isequal(s.type,'()')
   error('Invalid structure or cell subscript.');
end
for n=1:length(s.subs)
   if isa(s.subs{n},'extended')
      s.subs{n}=double(s.subs{n});
   end
end
if ~isa(o,'extended')
   o=extended(o);
end
if ~isa(v,'extended')
   v=extended(v);
end
d=size(o);
e=size(s.subs{1});
e=numel(s.subs) == 1 & isnumeric(s.subs{1}) & numel(e) == 2 & any(e == 1);
v=v.extended;
if e & numel(v) & numel(d) == 2 & d(1) < 2
   s.subs={':',1,s.subs{:}};
else
   s.subs={':',s.subs{:}};
end
switch numel(v)
case 0
   o.extended=subsasgn(o.extended,s,[]);
   if e
      [d e]=size(o.extended);
      o.extended=reshape(o.extended,[d 1 e]);
   end
   return
case 10
   o.extended=subsasgn(o.extended,s,0);
   if ~any(v(:))
      return
   end
   s.subs{1}=1;
   v=repmat(v,size(subsref(o.extended,s)));
   s.subs{1}=':';
end
o.extended=subsasgn(o.extended,s,v);
