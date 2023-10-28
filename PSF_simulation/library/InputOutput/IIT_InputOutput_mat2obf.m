function IIT_InputOutput_mat2obf(array)
% function IIT_InputOutput_mat2dat(array)
% Write the elements of the array into an .obf file (with the same name)
%
% array:        array to write
%

dim = ndims(array);

arrayname = inputname(1);
filename = strcat(arrayname,'.obf');

h = omas_bf_open(filename,1);
omas_bf_write(h,array);
omas_bf_close(h);

end
