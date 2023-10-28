function img = IIT_InputOutput_msr2mat(filename)
% IIT_InputOutput_msr2mat(filename)
% return the image contained into an .msr file into a verctor and save it
% into a .mat file (save name)
%
% filename:        name of the file to open
%

[h_image,info_file_image] = omas_bf_open(filename);
display(['Number of stacks stored in the file:' num2str(info_file_image.num_stacks)]);
img = omas_bf_read(h_image,1);
omas_bf_close(h_image);

outdataFilename = [filename(1:size(filename,2)-4) '.mat']
save(outdataFilename,'img');

end
