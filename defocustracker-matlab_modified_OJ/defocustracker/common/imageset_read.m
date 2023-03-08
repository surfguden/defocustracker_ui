function im = imageset_read(this,n,varargin)

im = give_image_raw(this,n);

if nargin==2
    if size(im,3)==3
        im = rgb2gray(im);
    end
    im = double(im);
end

if nargin>2
    if ~strcmp(varargin{1},'raw')
        disp('option not valid')
    end
end


function im = give_image_raw(this,n)
switch this.im_type
    case 'preloaded_imageset' %added functionality for preloaded images
        im=this.images{n};
    case 'DAVIS single'
        im0 = readimx(fullfile(this.path,this.images{n}));
        im = im0.Frames{1}.Components{1}.Planes{1}';
    case 'DAVIS double'
        n0 = round(n/2);
        fm = n-n0*2+2;
        im0 = readimx(fullfile(this.path,this.images{n0}));
        if fm==1
            im = im0.Frames{1}.Components{1}.Planes{1}';
        elseif fm==2
            im = im0.Frames{2}.Components{1}.Planes{1}';
        end
    case 'DAVIS'
        im0 = readimx(fullfile(this.path,this.images{1}));
        all_fr = length(im0.Frames);
        n_fr = mod(n,all_fr); n_fr = n_fr+double(n_fr == 0)*all_fr;
        n0 = (n-n_fr)/all_fr+1;
        im0 = readimx([this.path this.images{n0}]);
        im = im0.Frames{n_fr}.Components{1}.Planes{1}';
    case 'TIF'
        im = imread(fullfile(this.path,this.images{n}));
    case 'PNG'
        im = imread(fullfile(this.path,this.images{n}));
    case 'CZI'
        im = bfaOpenSingle(fullfile(this.path,this.images{1}),n);
end

if size(im,3)==3
    im = rgb2gray(im);
end
