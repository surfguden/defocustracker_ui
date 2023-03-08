function varargout = calibration_1_show(calib,varargin)

if nargin==1
    name_var = inputname(1);
elseif nargin==2
    name_var = varargin{1};
end

cal = calib.parameter;
cal.imresize = calib.training.imresize;

xy0 = (size(cal.images{1})+1)/2;                    % define origo
f = figure; set(f,'position',[100 50 900 600])  % initialize figure
set(f,'name',['model: ',name_var],'numbertitle','off')

% Plot: Signal-to-noise ratio vs. image number
axf2 = axes('Parent',f,'position',[.08 .61 .38 .36]);
plot(cal.z,cal.signal_to_noise,'.-');
hold on;
hpo = plot(cal.z(1),cal.signal_to_noise(1),'.','MarkerSize',15);     %Set the marker to the first image, which is at z=0
xlim([0 1]) 
xlabel('z/h'); ylabel('Signal-to-noise ratio');


% Plot: Similarity rate of change vs. image number
axf3 = axes('Parent',f,'position',[.6 .61 .38 .36]);
s_neighbor = diag(cal.similarity_map,-1)';
z_neighbor = cal.z(1:end-1)+diff(cal.z(1:2))/2;
plot(z_neighbor,s_neighbor,'.-');
hold on; 
hpo2 = plot(cal.z(1),s_neighbor(1),'.','MarkerSize',15);
xlim([0 1]) 
xlabel('z/h'); ylabel('Similarity');

% Plot: Correlation matrix surface
axes('Parent',f,'position',[.6 .1 .35 .35]);
if numel(cal.similarity_map)>1
    surf(cal.similarity_map)
else
    plot3(1,1,cal.similarity_map)
end
xlabel('j'), ylabel('i'), zlabel('Cm')
xlim([0, size(cal.similarity_map,2)+1]), ylim([0, size(cal.similarity_map,1)+1]), zlim([0 1])

% Plot: Particle image slider
axf = axes('Parent',f,'position',[.08 .16 .37 .34]);
% axes('Parent',f,'position',[.08 .61 .38 .36]);
% him = imagesc(cal.images{1}); 
him = imshow(cal.images{1},...
    [min(cal.images{1}(:)) max(cal.images{1}(:))+1],...
    'initialMagnification','fit'); 
hold on, daspect([1 1 1])
size_c = size(cal.images{1});
plot((size_c(2)+1)/2,(size_c(1)+1)/2,'+')
hpl = plot(cal.images_boundary{1}(1,:)*cal.imresize+xy0(2),cal.images_boundary{1}(2,:)*cal.imresize+xy0(1),'g');
set(gca,'YDir','normal',...
    'XTick',[1,(size_c(2)+1)/2,size_c(2)],...
    'YTick',[1,(size_c(1)+1)/2,size_c(1)])
xlabel('X'), ylabel('Y')
bl1 = uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.43 .015 .1 .04],'String',['1/',num2str(length(cal.images))],'Fontsize',11);
if length(cal.images)>1
    b = uicontrol('Parent',f,'Style','slider','units','normalized','SliderStep',[1 1]/(length(cal.images)-1),...
        'position',[.1 .02 .33 .04],'Value',1,'min',1,'max',length(cal.images));
    b.Callback = {@updateSlider,him,hpl,hpo,hpo2,xy0,cal.images_boundary,cal.images,bl1,axf,axf2,axf3,...
                    cal.imresize,cal.z,cal.signal_to_noise,z_neighbor,s_neighbor};
end

if nargout==1
    varargout{1} = s2;
end

function updateSlider(a,~,him,hpl,hpo,hpo2,xy0,image_boundary,im_c,bl1,axf,axf2,axf3,imr,cal_z,cal_snr,z_n,s_n)

n = round(a.Value);
if n==0, n=1; end
set(axf,'CLim',[min(im_c{n}(:)) max(im_c{n}(:)+1)])
set(him,'CData',im_c{n})
set(hpl,'XData',image_boundary{n}(1,:)*imr+xy0(2),...
    'YData',image_boundary{n}(2,:)*imr+xy0(1))
set(bl1,'String',[num2str(n),'/',num2str(length(im_c))])
% set(axf2)
cal_zh = [cal_z(:); cal_z(end)];
cal_snrh = [cal_snr(:); cal_snr(end)]; 
s_nh = [s_n(:); s_n(end)];
z_nh = [z_n(:); z_n(end)];

set(hpo,'XData',cal_zh(n),'YData',cal_snrh(n))
% set(axf3)
set(hpo2,'XData',z_nh(n),'YData',s_nh(n))

