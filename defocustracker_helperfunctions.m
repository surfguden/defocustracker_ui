classdef defocustracker_helperfunctions
   properties
      %Value {mustBeNumeric}
   end
   methods
       function median_imgs = median(obj,myfolder)
        %myfolder='C:\test\230131_temp data\mV1000\Repeat0'
        z_stack_imgs=dir(fullfile(myfolder, '**\*.tif'));
        folder={z_stack_imgs.folder};
        filenames={z_stack_imgs.name};
        fullpath=fullfile(folder, filenames);   
        %read all images
        imgs={};
        for i=1 : numel(fullpath)    
           imgs{i}=imread(fullpath{i});
        end
                       
        median_img=median(cat(3,imgs{:}),3);
        
        for i=1 : numel(fullpath)
            path=fullfile('C:\test\230131_temp data\median', filenames{i});
            %imtowrite=imgs{i}-median_img;
            median_imgs{i}=imgs{i}-median_img;
            %imwrite(imtowrite,path);
        end
       end
   end
end





