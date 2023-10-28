%Display array.
%
%    Marcel Leutenegger © 2.5.2005
%
function display(s)
if isequal(get(0,'FormatSpacing'),'compact')
   disp([inputname(1) ' =']);
else
   disp(' ');
   disp([inputname(1) ' =']);
   disp(' ');
end
if isempty(s)
   disp('     []');
   disp(' ');
else
   disp(double(s));
end
