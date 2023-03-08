function im_set = imageset_create(varargin)

im_type = [];
if nargin == 0
    [filename, pathname, ok] = uigetfile({'*.tif';'*.im7';...
        '*.png';'*.czi'},'Select set of images','MultiSelect','on');
    if ok>0
        if ~iscell(filename), im_name{1} = filename;
        else
            im_name = filename;
        end
        im_path = pathname;
    else
        im_name = [];
        im_path = [];
    end
elseif nargin == 1
    im_type = 'preloaded_imageset';
    im_name=varargin{1};
    im_path = [];
    n_frames = numel(varargin{1});
elseif nargin == 2
    im_path = varargin{1};
    im_name = varargin{2};
end
if ~isempty(im_name)
    if strcmp(im_name{1}(end-2:end),'tif')
        im_type = 'TIF';
        n_frames = length(im_name);
    elseif strcmp(im_name{1}(end-2:end),'im7')
        im0 = readimx(fullfile(im_path,im_name{1}));
        if length(im0.Frames)==1
            im_type = 'DAVIS single';
            n_frames = length(im_name);
        elseif  length(im0.Frames)==2
            im_type = 'DAVIS double';
            n_frames = length(im_name)*2;
        else
            im_type = 'DAVIS';
            n_frames = length(im_name)*length(im0.Frames);
        end
    elseif strcmp(im_name{1}(end-2:end),'czi')
        im_type = 'CZI';
        n_frames = bfaGetLength(fullfile(im_path,im_name{1}));
    elseif strcmp(im_name{1}(end-2:end),'png')
        im_type = 'PNG';
        n_frames = length(im_name);
    end
end
if isempty(im_type)
    disp('Empty image set or type not supported')
else
    im_set.path = im_path;
    im_set.images = im_name;
    im_set.n_frames = n_frames;
    im_set.im_type = im_type;
    im_set.type = 'imageset';
end
