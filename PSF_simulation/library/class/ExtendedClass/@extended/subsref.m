%Subscripted reference.
%
%    Marcel Leutenegger © 2.5.2005
%
function o=subsref(o,s)
if ~isequal(s.type,'()')
   error('Invalid structure or cell subscript.');
end
for n=1:length(s.subs)
   if isa(s.subs{n},'extended')
      s.subs{n}=double(s.subs{n});
   end
end
if isa(o,'extended')
   d=size(o.extended);
   e=size(s.subs{1});
   if length(d) == 3 & d(2) == 1 & d(3) > 1 & ...
      length(s.subs) == 1 & isnumeric(s.subs{1}) & ...
      length(e) == 2 & any(e == 1)
      s.subs={':',1,s.subs{:}};
   else
      s.subs={':',s.subs{:}};
   end
   o.extended=subsref(o.extended,s);
else
   o=subsref(o,s);
end
