function [dat, varargout] = process_0(model,imSet,selected_frames)

% This function has the purpose to identify particles in the images.
% This is done by differentiating the pixels based on a value threshold.
% It further defines the center coordinates of those detected particles.
% It also serves the purpose to allow tracking.

dat = dtracker_create('dataset');
RP = cell(1,imSet.n_frames);            % creates a cell array for the
% particle area boundaries for the
% respective number of frames

N_cal = length(selected_frames);


tic                                     % Starts the clock

for k = 1:N_cal                	% performs the particle detection and
    % creation of data for every image.
    
    IM1 = dtracker_show(imSet,selected_frames(k));
    [ellt, RP(selected_frames(k))] = single_frame_run(model,IM1,selected_frames(k));    % main function of the 2D-Tracking
    dat = [dat; ellt];	% stacks the 2D-tracked data from frame k to the existing dataset of the previously evaluated frames.
    tsec = toc; perc = (k/N_cal);   % defines the overall progress of the tracking
    disp([num2str(k),'/',num2str(N_cal),', ',secs2hms(tsec/k*(N_cal-k))])
end

end



function [ell, varargout] = single_frame_run(model,IM,n)
% process of detecting position of possible particles based on filtering
% threshold and choosen filters/preprocessing.

IM = preprocessing(model,IM);   	% The image gets filtered and preprocessed.

[X0,Y0,RP] = findsegm(model,IM);    % Finds the center and ellipses with
% pixel accuracy.

% subpixel accuracy - spatial correlation
%             if model.p.Depth_determination == 1 && ~isempty(model.cal.c.IC1)
%                 [X,Y,Z,Cm] = model.spatialcorr(IM,X0,Y0,RP);
%             elseif model.p.Depth_determination == 1 && isempty(model.cal.c.IC1)
%                 disp('Depth determination not possible: invalid calibration object')
%                 Z = X0*0;
%                 X = X0;
%                 Y = Y0;
%                 Cm = X0*0;
%             else
Z = X0*0;
X = X0;
Y = Y0;
Cm = X0*0;
%             end


t = ones(size(X))*n;

ell = dtracker_create('dataset',t,X,Y,Z,X*0,X*0,X*0,X*0,Cm);
% creates the dataset with      x,y,z,dx,dy,dz,fr,id,cm and
% scaling, type and metadata.
% The dataset from the respective frame includes the center coordinates of
% all detected ellipses in this frame.
% They need to be saved together with the frame number (variable t), to be
% identifiable later on.

% create RP
if nargout == 2
    RPnew{1}.B = RP.B;
    varargout{1} = RPnew;
end

end


function IM = preprocessing(model,IM)

if size(IM,3) == 3, IM = rgb2gray(IM); end  % Grayscaling, if needed.

% median filter:
if model.processing.median_filter > 1
    IM = medfilt2(IM,[model.processing.median_filter model.processing.median_filter]);
end
% gaussian filter:
if model.processing.gauss_filter > 1
    PSF = fspecial('gaussian',model.processing.gauss_filter,model.processing.gauss_filter);
    IM = imfilter(IM,PSF,'symmetric','conv');
end
end



function [X0,Y0,RP] = findsegm(model,IM)    % finds segments based on a brighter appearance in the image.

dumm = IM>model.processing.boundary_threshold; % Every value above the threshold gets a 1, the others a 0
dumm = imfill(dumm,'holes');    % fills areas surrounded by 1s up with 1s. This represents the part inside the brighter particle boundaries
[B,L] = bwboundaries(dumm,8);   % gives out the coordinates of every single boundary pixel
A = regionprops(L,'Area','Centroid','PixelIdxList');
AR = [A(:).Area]; ind = find(AR>model.processing.min_area);
X0 = ind*0;
Y0 = ind*0;
RP.L = L*0;
RP.B = {};
for k = 1:length(ind)
    X0(k) = A(ind(k)).Centroid(1);
    Y0(k) = A(ind(k)).Centroid(2);
    RP.L(A(ind(k)).PixelIdxList) = k;
    RP.B{k}(:,1) = B{ind(k)}(:,2);
    RP.B{k}(:,2) = B{ind(k)}(:,1);
end

end





%%

function time_string=secs2hms(time_in_secs)
time_string='';
nhours = 0;
nmins = 0;
if time_in_secs >= 3600
    nhours = floor(time_in_secs/3600);
    if nhours > 1
        hour_string = ' hours, ';
    else
        hour_string = ' hour, ';
    end
    time_string = [num2str(nhours) hour_string];
end
if time_in_secs >= 60
    nmins = floor((time_in_secs - 3600*nhours)/60);
    if nmins > 1
        minute_string = ' mins, ';
    else
        minute_string = ' min, ';
    end
    time_string = [time_string num2str(nmins) minute_string];
end
nsecs = time_in_secs - 3600*nhours - 60*nmins;
time_string = [time_string sprintf('%2.1f', nsecs) ' secs'];
end
