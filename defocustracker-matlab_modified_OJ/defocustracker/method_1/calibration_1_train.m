function cal_out = calibration_1_train(cal,imset,datset)

cal_out = cal;

iswin = ispc; % Check if system is windows

frames = datset.fr;
width = cal.training.imwidth;
height = cal.training.imheight;
x = datset.X;
y = datset.Y;
stack.z = datset.Z;
stack.images = cell(1,length(frames));
for ii = 1:length(frames)
    im = imageset_read(imset,frames(ii));
    imc = imageset_crop(im,x(ii),y(ii),[width, height]);
    stack.images{ii} = imc;
end


parameters = cal.training;

median_filter = parameters.median_filter;
gauss_filter = parameters.gauss_filter;
smoothing = parameters.smoothing;
n_interp_images = parameters.n_interp_images;
boundary_threshold = parameters.boundary_threshold;

n_images = length(stack.images);
stack_index = 1:n_images;

disp('Training started ...')
tic

% Internal parameters for calibration
pad_im_t = 2;

%% Find im_c and im_t and signal-to-noise ratio
background = stack_index*0;
noise = stack_index*0;
im_c = cell(1,n_images);
im_t = im_c;
n_pixels = numel(stack.images{1});
im_ca = zeros(n_pixels,n_images);
for ii = 1:n_images
    im = stack.images{ii};
    if median_filter>2
        im =  medfilt2(im,[median_filter median_filter],'symmetric');
    end
    if gauss_filter>1
        PSF = fspecial('gaussian',gauss_filter,gauss_filter);
        im = imfilter(im,PSF,'symmetric','conv');
    end
    im_c{ii} = im;
    im_t{ii} = padarray(im,pad_im_t,'symmetric');
    background(ii) = median([im(1,:), im(end,:), im(:,1)', im(:,end)']);
    noise(ii) = std([im_c{ii}(1,:), im_c{ii}(end,:), im_c{ii}(:,1)', im_c{ii}(:,end)']);
    im_ca(:,ii) = im_c{ii}(:);
end
background = median(background);
noise = median(noise);



%% Smoothing and interpolation

if ~(smoothing==0 && n_interp_images==0)
    
    [n_I, n_J] = size(im_c{1});
    n_pixels = numel(im_c{1});
    n_stack = length(im_c);
    ind_stack = stack.z;
    
    im_new = zeros(n_pixels,n_stack);
    for ii = 1:n_stack
        im_new(:,ii) = im_c{ii}(:);
    end
    if smoothing~=0
        im_orig = im_new;
        im_new = zeros(n_pixels,n_stack);
        for ii = 1:n_pixels
            im_new(ii,:) = smooth(ind_stack,im_orig(ii,:),smoothing,'loess');
        end
    end
    if n_interp_images~=0
        im_orig = im_new;
        im_new = zeros(n_pixels,n_interp_images);
        ind_new = linspace(0,1,n_interp_images);
        for ii = 1:n_pixels
            im_new(ii,:) = interp1(ind_stack,im_orig(ii,:),ind_new);           
        end
    end
    im0 = zeros(n_I, n_J);
    im_c = cell(1,size(im_new,2));
    im_t = cell(1,size(im_new,2));
    for ii = 1:length(im_c)
        im0(:) = im_new(:,ii);
        im_c{ii} = im0;
        im_t{ii} = padarray(im0,pad_im_t,'symmetric');
    end
    
    n_images = length(im_c);
end


%% reduce images
if parameters.imresize~=1
    for ii = 1:length(im_c)
        im_c{ii} = imresize(im_c{ii},parameters.imresize);
        im_t{ii} = padarray(im_c{ii},pad_im_t,'symmetric');
    end
end

[n_i, n_j] = size(im_c{ii});
if mod(n_i,2)==0
    for ii = 1:length(im_c)
        im_c{ii} =  (im_c{ii}(2:end,:)+im_c{ii}(1:end-1,:))/2;
        im_t{ii} = padarray(im_c{ii},pad_im_t,'symmetric');
    end
end
if mod(n_j,2)==0
    for ii = 1:length(im_c)
        im_c{ii} =  (im_c{ii}(:,2:end)+im_c{ii}(:,1:end-1))/2;
        im_t{ii} = padarray(im_c{ii},pad_im_t,'symmetric');
    end
end


%% Find particle image borders
image_boundary = cell(size(im_c));
image_area = cell(size(im_c));
boundary_threshold = max([noise 1])*boundary_threshold;
sizeim_c = size(im_c{1});
signal = zeros(size(im_c));

% lrbt_im_c = zeros(4,length(im_c));
for ii = 1:n_images
    [~,im_segm] = bwboundaries(abs(im_c{ii}-background)>boundary_threshold,8);
    segm_data = regionprops(im_segm,'PixelIdxList');
    id_x = []; id_y = [];
    if ~isempty(segm_data)
        for jj = 1:length(segm_data)
            if length(segm_data(jj).PixelIdxList)>10
                [id_y_dum, id_x_dum] = ind2sub(size(im_segm),segm_data(jj).PixelIdxList);
                id_x = [id_x; id_x_dum]; id_y = [id_y; id_y_dum];
            end
        end
    end
    
    if ~isempty(id_x)
        ind_convex_hull = boundary(id_x,id_y,0);
        image_boundary{ii}(1,:) = id_x(ind_convex_hull)-(sizeim_c(2)+1)/2;
        image_boundary{ii}(2,:) = id_y(ind_convex_hull)-(sizeim_c(1)+1)/2;
        image_area{ii}(1,:) = id_x-(sizeim_c(2)+1)/2;
        image_area{ii}(2,:) = id_y-(sizeim_c(1)+1)/2;
        linearimnd = sub2ind(sizeim_c,id_y,id_x);
        signal(ii) = mean(im_c{ii}(linearimnd));
    else
        image_boundary{ii}(1,:) = [-.5 -.5 .5 .5 -.5]*(size(im_c{ii},2)-1);
        image_boundary{ii}(2,:) = [-.5 .5 .5 -.5 -.5]*(size(im_c{ii},1)-1);
        [dumx, dumy] = meshgrid(image_boundary{ii}(1,:),image_boundary{ii}(2,:));
        image_area{ii}(1,:) = dumx(:)';
        image_area{ii}(2,:) = dumy(:)';
        signal(ii) = mean(im_c{ii}(:));
    end
end

signal_to_noise = (signal-background)/noise;

%% Find cal_05: pixel with a signal through the whole scan
global_area = zeros(size(im_c{1}));
idx = []; idy = [];

for ii = 1:length(image_area)
    idx = [idx, image_area{ii}(1,:)];
    idy = [idy, image_area{ii}(2,:)];
end
idx = idx+(size(global_area,2)+1)/2;
idy = idy+(size(global_area,1)+1)/2;
linearimnd = sub2ind(size(global_area),idy,idx);
linearimnd = unique(linearimnd);
global_area(linearimnd) = 1;



%% Create correlation matrix
similarity_map = zeros(length(im_c),length(im_t));
for ii = 1:length(im_c)
    for jj = 1:length(im_t)
        dum = normxcorr2_general(im_c{ii},im_t{jj},iswin);
        similarity_map(ii,jj) = max(dum(:));
    end
end


%% Update similarity map

if length(im_c)==length(im_t)
    
    s_new = similarity_map(1,:)*0;
    for ii = 1:size(similarity_map,2)
        indii = [-2 -1 1 2];
        switch ii
            case 1, indii = [1 2 3];
            case 2, indii = [-1 1 2];
            case size(similarity_map,2)-1, indii = [-2 -1 1];
            case size(similarity_map,2), indii = [-3 -2 -1];
        end
        x0 = indii+ii;
        y0 = similarity_map(ii,x0);
        
        system_matrix = [x0'.^2-2*ii*x0' x0'.*0+1];
        coefficients = (system_matrix'*system_matrix)\(system_matrix'*y0');
        coeffc = [coefficients(1) -2*coefficients(1)*ii coefficients(2)];
        s_new(ii) = polyval(coeffc,ii);
    end
    id_ij = logical(eye(size(similarity_map)));
    similarity_map(id_ij) = s_new;
else
    similarity_old = similarity_map;
    n_map = size(similarity_old,1);
    ind_old = 1:size(similarity_old,2);
    ind_map = linspace(1,size(similarity_old,2),n_map);
    similarity_map = zeros(n_map);
    for ii = 1:n_map
        similarity_map(ii,:) = interp1(ind_old,...
            similarity_old(ii,:),ind_map,'spline');
    end
end

%%
if parameters.imresize~=1
    for ii = 1:length(im_c)
        image_boundary{ii} = image_boundary{ii}/parameters.imresize;
    end
end

mean_area = 0; 
for ii = 1:length(image_area)
    mean_area = mean_area + size(image_area{ii},2); 
end
mean_area = mean_area/length(image_area);


%% cal struct: data used in the frame run

% calibration data
cal_out.parameter.images = im_c;
cal_out.parameter.images_boundary = image_boundary;
cal_out.parameter.images_area = image_area;
cal_out.parameter.z = linspace(0,1,length(im_c));
cal_out.parameter.signal_to_noise = signal_to_noise;
cal_out.parameter.similarity_map = similarity_map;
cal_out.parameter.global_area = global_area;
cal_out.parameter.mean_area = mean_area;
cal_out.parameter.n_cal = length(im_c);


% Display total time
tot_time = toc;
timestr = secs2hms(tot_time);
disp(['Training done! Total Time: ',timestr])


