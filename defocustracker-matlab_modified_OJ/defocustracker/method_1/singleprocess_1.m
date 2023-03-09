function [dat, eval_time, im] = singleprocess_1(im,cal,par)


tic

% Check if system is windows
iswin = ispc;

% 2D GDPT to be implemented (see code in the end)
%%
% Applying image filter
if cal.median_filter>2
    im =  double(medfilt2(uint16(im),[cal.median_filter cal.median_filter],'symmetric'));
end
if cal.gauss_filter>1
    PSF = fspecial('gaussian',cal.gauss_filter,cal.gauss_filter);
    im = imfilter(im,PSF,'symmetric','conv');
end
if cal.imresize~=1
    im = imresize(im,cal.imresize);
end
%%

% The process is repeated to find more particle. It works only if the
% walking step is on
if par.walking_step==0, par.n_iteration=1; end
% Standard evaluation require 1 iteration (par.n_iteration = 1);
x3 = []; y3 = []; z3 = []; c3 = [];

for nn = 1:par.n_iteration
    [xr, yr, zr, cr, im] = find_positions(im,par,cal,iswin);
    
    x3 = [x3, xr]; y3 = [y3, yr]; z3 = [z3, zr]; c3 = [c3, cr];
end

% Without sub-image interpolation the output is not normalized!
% ONLY ADVAMCED USERS SHOULD TURN OFF SUBIMAGE INTERPOLATION
if par.subimage_interpolation==1
   flag = z3~=1 & z3~=cal.n_cal; 
   x3 = x3(flag); y3 = y3(flag); z3 = z3(flag); c3 = c3(flag);
   z3 = (z3-1)/(cal.n_cal-1);
end

% Rescale x and y values in case of image resize
if cal.imresize~=1
    x3 = x3/cal.imresize;
    y3 = y3/cal.imresize;
end

% Write final DataSet
dat = dataset_create(x3*0,x3,y3,z3,x3*0,x3*0,x3*0,x3*0,c3);
eval_time = toc;

function [x3, y3, z3, c3, im] = find_positions(im,par,cal,iswin)


% initialize parameters
pad_im_t = 4; % this must be the same also in the calibration
size_t = size(cal.images{1})+2*pad_im_t; % size of target image
size_c = size(cal.images{1}); % size of calibration image
stack_length = length(cal.images); % total number of stack images

%% STEP 1 - FIND PARTICLE POSITIONS - 
% guess xy-position and calculate global similarity maps from a
% reduced set of stack images
x1 = []; y1 = []; c1 = [];
if canUseGPU()==1
    im=gpuArray(im);
end

% find index of reduced number of stack images used in the first step
if length(par.images_in_guess_step)==1 
    ind_stack = round(linspace(1,stack_length,par.images_in_guess_step*2+1));
    ind_stack = ind_stack(2:2:end-1);
else
    ind_stack = par.images_in_guess_step;
    disp('selcted custom indstack')
end

ccn = zeros(size(im,1),size(im,2),length(ind_stack));
for ii = 1:length(ind_stack)
    cc = normxcorr2(cal.images{ind_stack(ii)},im);
    % cut correlation to the size of image
    cc_im = size(cc)-size(im);
    ccn(:,:,ii) = cc(1+ceil(cc_im(1)/2):end-floor(cc_im(1)/2),1+ceil(cc_im(2)/2):end-floor(cc_im(2)/2));
    % find peaks, i.e. candidate particles
    [xn,yn,cn] = find_peaks_2d(ccn(:,:,ii),par.cm_guess);
    x1 = [x1; xn]; y1 = [y1; yn];  c1 = [c1; cn];
end

% remove particles too close to the border
flag_out = x1<2+size_t(2)/2 | x1>size(im,2)-size_t(2)/2-1 |...
    y1<2+size_t(1)/2 | y1>size(im,1)-size_t(1)/2-1;
x1(flag_out) = []; y1(flag_out) = []; c1(flag_out) = [];
if canUseGPU()==1
    im=gather(im);
end

%% STEP 2 - GUESSING STEP - find a first guess of z position from the reduced
x2 = []; y2 = []; z2 = []; c2 = []; ts2 = [];

% indices of target image
ii_t = round(-(size_t(1))/2)+1:round((size_t(1))/2)-1; 
jj_t = round(-(size_t(2))/2)+1:round((size_t(2))/2)-1;
% circle of no overlapping zone
circx = par.no_overlap_radius*cos(linspace(0,2*pi,16));
circy = par.no_overlap_radius*sin(linspace(0,2*pi,16));

% iteration over all the detected points in the FIRST STEP
sim_map0 = cal.similarity_map(ind_stack,:); % similarity map reduced
while ~isempty(c1)
    nn = find(c1==max(c1),1);
    
    sim_vect = squeeze(ccn(y1(nn),x1(nn),:));
    [new_ind, top_sim] = guess_stack(sim_map0,sim_vect);
    ind_vect = ind_stack;
    % sim_vect: vector of similarity values
    % ind_vect: index in the stack of the values of sim_vect
    % new_ind: current index of best match
    
    im_t = im(ii_t+round(y1(nn)),jj_t+round(x1(nn))); % target image
    % --- the guess is continued until new_ind is contained in ind_vect
    is_match_vect = 0;
    while ~any(ind_vect == new_ind)
        ind_vect = [ind_vect, new_ind];
        dum = normxcorr2_general(cal.images{new_ind},im_t,iswin);
        sim_vect = [sim_vect; max(dum(:))];
        sim_map = cal.similarity_map(ind_vect,:);
        [new_ind, top_sim, match_vect] = guess_stack(sim_map,sim_vect);
        is_match_vect = 1;
    end
    ts2 = [ts2, top_sim]; % IS THIS USEFUL??? Should be exchanged with c2?
    
    if par.walking_step==0 && par.subimage_interpolation==1 && is_match_vect==1
        if new_ind>1 && new_ind<length(match_vect)
            sim_vect2 = match_vect(new_ind-1:new_ind+1);
            new_ind = new_ind+(sim_vect2(1)-sim_vect2(3))/...
                (2*sim_vect2(1)-4*sim_vect2(2)+2*sim_vect2(3));
        end
    end

    x2 = [x2, x1(nn)]; y2 = [y2, y1(nn)]; z2 = [z2, new_ind]; c2 = [c2, c1(nn)];
    % flag out particles inside the no-overlapping circle
    flag_out = inpolygon(x1,y1,x1(nn)+circx,y1(nn)+circy);
    x1(flag_out) = []; y1(flag_out) = []; c1(flag_out) = [];
end

if par.walking_step==0 
    flag_out = isnan(x2) | isnan(y2) | isnan(z2) | c2<par.cm_guess;
    x3 = x2(~flag_out);
    y3 = y2(~flag_out);
    z3 = z2(~flag_out);
    c3 = c2(~flag_out);
    return
end


%%% CHECK!!!
% flag_in = ts2>.8;
% x2 = x2(flag_in); y2 = y2(flag_in); z2 = z2(flag_in); c2 = c2(flag_in);
%% STEP 3 - WALKING STEP - find final xyz position with walking!
x3 = NaN(size(x2)); y3 = x3; z3 = x3; c3 = x3; ts3 = ts2;
nxcorrs = c2*0; % debugging parameter to count number of correlation

% DEBUGGING: TRICK-OVERLAP are used to try to improve the the detection of
% overlapping particles. The code belongin to one specific trick is marked
% with TRICK-OVERLAP and the respective number

%%% TRICK-OVERLAP-1 -- masking pixel outside the particle.
%%% Drawback: does not work properly is the particle is not well centered
%%% after the SECOND STEP
% this block create a mask for the target image. The mask indicate pixel in
% which there is some signal throught the stack.
im_t_bk = padarray(cal.global_area,[pad_im_t pad_im_t]);
im_t_bkout = double(im_t_bk==0);
%%%%%%%%

%%% TRICK-OVERLAP-2 -- re-run the SECOND STEP if
im_ghost = im*0;
min_ghost = 10;
no_overlap = 0;
%%%%%%%%

%%% TRICK-OVERLAP-3 -- "Delete" the particle from the image,
% see end of for loop


for ii = 1:length(x2)

    xx = x2(ii); yy = y2(ii); new_ind = z2(ii); cs = c2(ii);

    nxcorr = 0;
    % find target image
    im_t = im(ii_t+round(yy),jj_t+round(xx));
    %%% TRICK-OVERLAP-1
    %     bkg2 = median([im_t(1,:),im_t(end,:),im_t(:,1)',im_t(:,end)']);
    bkg2 = median(im_t(:));
    im_t = im_t.*im_t_bk+im_t_bkout*bkg2;
    
    %%%
    
    %%% TRICK-OVERLAP-2
    im_t_ghost =  im_ghost(ii_t+round(yy),jj_t+round(xx));
    if sum(im_t_ghost(:))>min_ghost && cs<par.cm_final
        sim_vect = ind_stack'*0;
        for nn = 1:length(ind_stack)
            cc = normxcorr2(cal.images{ind_stack(nn)},im_t);
            sim_vect(nn) = max(cc(:));
            sim_map = cal.similarity_map(ind_stack,:);
            [new_ind, top_sim] = guess_stack(sim_map,sim_vect);
        end
        ts3(ii) = top_sim;
        no_overlap = no_overlap+1;
    end
    %%%
    
    field = normxcorr2(cal.images{new_ind},im_t); nxcorr = nxcorr+1;
    
    % jump to next iteration if no proper xcorr
    if max(field(:))<par.cm_guess, continue, end
    
    [iy, ix] = find(field==max(field(:)),1);
    xx = xx - (size_t(2)+size_c(2))/2  + ix + (log(field(iy,ix-1))-log(field(iy,ix+1))) /...
        ( 2*log(field(iy,ix-1)) - 4* log(field(iy,ix)) + 2* log(field(iy,ix+1)));
    yy = yy - (size_t(1)+size_c(1))/2 + iy + (log(field(iy-1,ix))-log(field(iy+1,ix))) /...
        ( 2*log(field(iy-1,ix)) - 4* log(field(iy,ix)) + 2* log(field(iy+1,ix)));
    
    % find a new target image with subpixel accuracy based on xx and yy
    if imag(xx)~=0 || imag(yy)~=0, continue, end
    im_t = mcrop(im,xx,yy,size_t);
    if isempty(im_t), continue, end
    
    %%% TRICK-OVERLAP-1
    im_t = im_t.*im_t_bk+im_t_bkout*bkg2; %%%
    
    dum = normxcorr2_general(cal.images{new_ind},im_t,iswin); nxcorr = nxcorr+1;
    maxdum = max(dum(:)); %  maxdum = dum(CimJ(1),CimJ(2));
    
    % --- Start "walking" to detect best maximum
    % initialize walking
    if new_ind==1
        ind_vect2 = [new_ind, new_ind+1, new_ind+2];
        dum1 = normxcorr2_general(cal.images{ind_vect2(2)},im_t,iswin); nxcorr = nxcorr+1;
        dum2 = normxcorr2_general(cal.images{ind_vect2(3)},im_t,iswin); nxcorr = nxcorr+1;
        sim_vect2 = [maxdum, max(dum1(:)), max(dum2(:))];
        %         sim_vect2 = [maxdum, dum1(CimJ(1),CimJ(2)), dum2(CimJ(1),CimJ(2))];
        dum123 = {dum,dum1,dum2};
    elseif new_ind==stack_length
        ind_vect2 = [new_ind-2, new_ind-1, new_ind];
        dum1 = normxcorr2_general(cal.images{ind_vect2(1)},im_t,iswin); nxcorr = nxcorr+1;
        dum2 = normxcorr2_general(cal.images{ind_vect2(2)},im_t,iswin); nxcorr = nxcorr+1;  
        sim_vect2 = [max(dum1(:)), max(dum2(:)), maxdum];
        %         sim_vect2 = [dum1(CimJ(1),CimJ(2)), dum2(CimJ(1),CimJ(2)), maxdum];
        dum123 = {dum1,dum2,dum};
    else
        ind_vect2 = [new_ind-1, new_ind, new_ind+1];
        dum1 = normxcorr2_general(cal.images{ind_vect2(1)},im_t,iswin); nxcorr = nxcorr+1;
        dum2 = normxcorr2_general(cal.images{ind_vect2(3)},im_t,iswin); nxcorr = nxcorr+1;         
        sim_vect2 = [max(dum1(:)), maxdum, max(dum2(:))];
        %         sim_vect2 = [dum1(CimJ(1),CimJ(2)), maxdum, dum2(CimJ(1),CimJ(2))];
        dum123 = {dum1,dum,dum2};
    end
    [maxdum, new_ind2] = max(sim_vect2);
    % (???) rename dum to put the final dum as the good one
    addim_ter = 2; whileBreak = false;
    
    % start walking
    while new_ind2~=2
        if new_ind2==1 && ind_vect2(1)>1
            ind_vect2 = [ind_vect2(1)-1, ind_vect2(1:2)];
            dum = normxcorr2_general(cal.images{ind_vect2(1)},im_t,iswin); nxcorr = nxcorr+1;   
            sim_vect2 = [max(dum(:)), sim_vect2(1:2)];
            %             sim_vect2 = [dum(CimJ(1),CimJ(2)), sim_vect2(1:2)];
            dum123 = {dum, dum123{1:2}};
            addim_ter = addim_ter+1;
        elseif new_ind2==3 && ind_vect2(3)<stack_length
            ind_vect2 = [ind_vect2(2:3), ind_vect2(3)+1];
            dum = normxcorr2_general(cal.images{ind_vect2(3)},im_t,iswin); nxcorr = nxcorr+1;          
            sim_vect2 = [sim_vect2(2:3), max(dum(:))];
            %             sim_vect2 = [sim_vect2(2:3), dum(CimJ(1),CimJ(2))];
            dum123 = {dum123{2:3}, dum};
            addim_ter = addim_ter+1;
        else
            whileBreak = true;
            [maxdum, new_ind2] = max(sim_vect2);
            break
        end
        [maxdum, new_ind2] = max(sim_vect2);
    end
    
    % --- Parabolic three-point fit to get final z with sub-image accuracy
    % COMMENT TO HAVE THE TEST WITHOUT sub-image accuracy
    % HINT: consider a 3rd polynomial fit?
    if ~whileBreak && par.subimage_interpolation == 1
        new_ind = ind_vect2(2);
        zz = new_ind+(sim_vect2(1)-sim_vect2(3))/...
            (2*sim_vect2(1)-4*sim_vect2(2)+2*sim_vect2(3));
        field = dum123{2};
    else
        new_ind = ind_vect2(new_ind2);
        field = dum123{new_ind2};
        zz = new_ind;
    end
    size_c = size(cal.images{new_ind});
    
    % Gaussian 3-point fit for final x-y with sub-pixel accuracy
    [iy, ix] = find(field==max(field(:)));
    if length(iy)==1
        xxfin = xx-(size_t(2)+size_c(2))/2  + ix + (log(field(iy,ix-1))-log(field(iy,ix+1))) /...
            ( 2*log(field(iy,ix-1)) - 4* log(field(iy,ix)) + 2* log(field(iy,ix+1)));
        yyfin = yy-(size_t(1)+size_c(1))/2 + iy + (log(field(iy-1,ix))-log(field(iy+1,ix))) /...
            ( 2*log(field(iy-1,ix)) - 4* log(field(iy,ix)) + 2* log(field(iy+1,ix)));
        if isnan(xxfin) || isnan(yyfin)
            %            disp('vaffanculo!') % hahahahahahahaha, mi piace
            xxfin = xx-(size_t(2)+size_c(2))/2 + sum(ix)/2;
            yyfin = yy-(size_t(1)+size_c(1))/2 + sum(iy)/2;
            %             field = normxcorr2(cal.images{new_ind},im_t);
            %             [iy, ix] = find(field==max(field(:)));
            %             xxfin = xx-(size_t(2)+size_c(2))/2  + ix + (log(field(iy,ix-1))-log(field(iy,ix+1))) /...
            %                 ( 2*log(field(iy,ix-1)) - 4* log(field(iy,ix)) + 2* log(field(iy,ix+1)));
            %             yyfin = yy-(size_t(1)+size_c(1))/2 + iy + (log(field(iy-1,ix))-log(field(iy+1,ix))) /...
            %                 ( 2*log(field(iy-1,ix)) - 4* log(field(iy,ix)) + 2* log(field(iy+1,ix)));
        end
    elseif length(iy)==2
        xxfin = xx-(size_t(2)+size_c(2))/2 + sum(ix)/2;
        yyfin = yy-(size_t(1)+size_c(1))/2 + sum(iy)/2;
    else
        xxfin = xx-(size_t(2)+size_c(2))/2 + ix(1);
        yyfin = yy-(size_t(1)+size_c(1))/2 + iy(1);
    end
    
    %%%% TRICK-OVERLAP-3!!!
    if ~isnan(zz) && maxdum>par.cm_final
        ii2 = round(yy+cal.images_area{round(zz)}(2,:));
        jj2 = round(xx+cal.images_area{round(zz)}(1,:));
        linearimnd = sub2ind(size(im),ii2,jj2);
        im(linearimnd) = bkg2;
        im_ghost(linearimnd) = 1;
    end
    %%%
    
    % create the final output of STEP THREE
    
    x3(ii) = xxfin; y3(ii) = yyfin;  z3(ii) = zz;  c3(ii) = maxdum;
    nxcorrs(ii) = nxcorr;
        
    %     catch
    %         keyboard
    %         continue
    %     end
end

%%

% mask out nan
flag_out = isnan(x3) | isnan(y3) | isnan(z3) | c3<par.cm_final;
x3(flag_out) = []; y3(flag_out) = []; z3(flag_out) = [];
c3(flag_out) = []; nxcorrs(flag_out) = []; ts3(flag_out) = [];



% disp(['number overlapping particles: ',num2str(no_overlap)])

% figure(gcf), clf, imagesc(im), hold on, plot(x2,y2,'+y',x3,y3,'+r')


function [new_ind, top_sim, match_vect] = guess_stack(sim_map,sim_vect)

% % Least square method (not robust against intensity variations)
% match_vect = sqrt(sum((sim_map-repmat(sim_vect,1,size(sim_map,2))).^2))/length(sim_vect);
% [~, new_ind] = min(match_vect);


% Normalized cross correlation method
match_vect = zeros(size(sim_map(1,:)));
mean_sim = mean(sim_vect);
for jj = 1:length(match_vect)
    sim_vect_jj = sim_map(:,jj);
    mean_sim_jj = mean(sim_vect_jj);
    
    % Discrete normalized 1-D cross correlation
    match_vect(jj) = sum((sim_vect_jj-mean_sim_jj).*(sim_vect-mean_sim))/...
        sqrt(sum((sim_vect_jj-mean_sim_jj).^2)*sum((sim_vect-mean_sim).^2));
end
[top_sim, new_ind] = max(match_vect);
% plot(1:length(match_vect),match_vect,'.-',new_ind,match_vect(new_ind),'*')


function imc = mcrop(im0,xo,yo,size_t)

xj = (-(size_t(2)-1)/2:(size_t(2)-1)/2)+xo;
yi = (-(size_t(1)-1)/2:(size_t(1)-1)/2)+yo;

dd = zeros(1,4); xx = dd; yy = dd;
IMs = cell(1,4);


if floor(xj(1))<1 || floor(yi(1))<1 ||...
        ceil(xj(end))>size(im0,2) || ceil(yi(end))>size(im0,1)
    imc = []; return
end


IMs{1} = im0(floor(yi),floor(xj));
xx(1) = xj(1)-floor(xj(1)); yy(1) = yi(1)-floor(yi(1));
IMs{2} = im0(ceil(yi),floor(xj));
xx(2) = xj(1)-floor(xj(1)); yy(2) = yi(1)-ceil(yi(1));
IMs{3} = im0(floor(yi),ceil(xj));
xx(3) = xj(1)-ceil(xj(1)); yy(3) = yi(1)-floor(yi(1));
IMs{4} = im0(ceil(yi),ceil(xj));
xx(4) = xj(1)-ceil(xj(1)); yy(4) = yi(1)-ceil(yi(1));
dd = sqrt(xx.^2+yy.^2);
[~, ii] = sort(dd);
xx = abs(xx(ii(1:3))); yy = abs(yy(ii(1:3)));
a = -xx(1)-yy(1)+1; b = xx(1); c = yy(1);
imc = IMs{ii(1)}*a+IMs{ii(2)}*b+IMs{ii(3)}*c;



%%%%
% GDPT for calibration stack with 1 image -- 2D PTV (see code in the end)
% to be implemented

% if length(cal.nArray)==1
%     imM = double(medfilt2(uint16(im),[cal.median_filter cal.median_filter],'symmetric'));
%     im_ts = imresize(imM,ares);
%     cc = normxcorr2_mex(cal.images{1},im_ts);
%     imJi = size(cc)-size(im_ts);
%     field = cc(1+ceil(imJi(1)/2):end-floor(imJi(1)/2),1+ceil(imJi(2)/2):end-floor(imJi(2)/2));
%     [x1, y1, Cm] = find_peaks_2d(field,par.cm_guess);
%     x = x1*0; y = y1*0;
%     for n = 1:length(x1)
%         ix = round(x1(n));
%         iy = round(y1(n));
%         x(n) = ix + (log(field(iy,ix-1))-log(field(iy,ix+1))) /...
%             ( 2*log(field(iy,ix-1)) - 4* log(field(iy,ix)) + 2* log(field(iy,ix+1)));
%         y(n) = iy + (log(field(iy-1,ix))-log(field(iy+1,ix))) /...
%             ( 2*log(field(iy-1,ix)) - 4* log(field(iy,ix)) + 2* log(field(iy+1,ix)));
%     end
%     flag_out = isnan(x);
%     x(flag_out) = []; y(flag_out) = []; Cm(flag_out) = [];
%     dat = cpoints(x,y,x*0+1,x*0+n,x*0,x*0,x*0,Cm,Cm*0);
%     eval_time = toc;
%     return
% end
