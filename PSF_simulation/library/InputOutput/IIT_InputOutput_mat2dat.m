function IIT_InputOutput_mat2dat(array)
% function IIT_InputOutput_mat2dat(array)
% Write the elements of array A to a binary file (with the same name)
% in column order and with the same precision of the array
%
% array:        array to write
%

arrayname = inputname(1);
prec = class(array)
filename = strcat(arrayname,'.dat');

fid = fopen(filename, 'w');
fwrite(fid, array, prec);
fclose(fid);

end

