function IIT_InputOutput_tiff2dat(filename)
% function [output_args] = IIT_InputOutput_tiff2dat(filename)
% Converts a tiff format file into a raw data format file.
% If filename is a directory converts all the file in the directory.

if isdir(filename)
    % read directory
    list = dir(filename);
    n_file = size(list,1);
    for i = 3:n_file
        disp(['Converting ', list(i).name]);
        
        try
            data2D = imread([filename,'\',list(i).name]);
            dot = strfind(list(i).name,'.');
            last_dot = dot(length(dot));
            format = list(i).name(last_dot+1:end);
            name = list(i).name(1:last_dot-1);
            
            filenameout = strcat(name, '.dat');
            
            % write file
            fid = fopen([filename,'\',filenameout], 'w');
            fwrite(fid, data2D, 'uint16', 0, 'l'); 
            fclose(fid);
        catch err
             disp(err.message);
        end
        
        
    end
else
    % read single file
    disp(['Converting  ', filename]);
    
    try
        data2D = imread(filename);
        dot = strfind(filename,'.');
        last_dot = dot(length(dot));
        format = filename(last_dot+1:end);
        name = filename(1:last_dot-1);
        
        filenameout = strcat(name, '.dat');
        
        % write file
        fid = fopen(filenameout, 'w');
        fwrite(fid, data2D, 'uint16', 0, 'l');
        fclose(fid);
    catch err
        disp(err.message);
    end
    
    
end
end

