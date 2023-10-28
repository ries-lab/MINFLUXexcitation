
%converts any .nd2 file into Imspector .obf file only when the structure is
%multible stacks within one window. No z-series!

%function = conv_nd2_msr('filename')
close all,
clear all,

directory='C:\Documents and Settings\bharke\My Documents\Bharke\Projects\Nikon\Experiments'; %directory, where the nd2 files and where the obf will be placed

files=dir(strcat(directory,'\*.nd2'));
for k=1:size(files,1)
    %file='01' %filename of the nd2 fie without ending
    %ending='.nd2'
    %filename=strcat(file,ending);
    filename=files(k).name;
    data=bfopen(strcat(directory,'\',filename));
    %data=bfopen('01.nd2');
    series1 = data{1, 1}; %for the 1st series
    series1_numPlanes = size(series1, 1); %gives number of planes of the 1st series
    pixel_size=str2num(data{1,2}.get('dCalibration')); %gets pixel size in nm from the metadata
    %series1
    data_01=cell(1,series1_numPlanes);
    for i=1:series1_numPlanes
        %data_01_data=sprintf('series1_stack%01d',i);
        data_01{1,i}=series1{i,1};
        %data_01_label=sprintf('series1_label%01d',i);
        %data_01_label(i)=series1{i,2};

        figure(k)
        subplot(1,series1_numPlanes,i)
        x_axes=0:pixel_size:size(data_01{1,i},1)*pixel_size; %giving pixel physical coordinates
        y_axes=x_axes; %only for quadratic pixel sizes
        imagesc(x_axes,y_axes,data_01{1,i});
        %title(data_01_label(i));
        colormap hot
        hold on
        axis image
    end
hold off


    pixel_num=size(data_01{1,1},1);
    %building the info struct for Imspector. It contains physical length of
    %the stack, stack names, class, etc... at the ,moment, the infos are
    %for all stacks within a measurement the same
    info=cell(1,series1_numPlanes);
    for j=1:series1_numPlanes
 
        info{1,j}=struct('name','Nikon','len',[pixel_size*pixel_num pixel_size*pixel_num],'cls','double'); %the metainfo for all stacks is the same
        
    end
    imspector_filename=strcat(directory,'\',filename,'.obf');
    omas_bf_write(imspector_filename,data_01,info);
    clear info;
    clear data;
    clear data_01;
    
end

