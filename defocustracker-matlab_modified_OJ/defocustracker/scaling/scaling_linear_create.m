function varargout = scaling_linear_create(varargin)

if nargout==1 && nargin==0
    Scaling.type = 'linear';
    Scaling.X_to_x = 1;
    Scaling.Y_to_y = 1;
    Scaling.Z_to_z = 1;
    Scaling.dt = 1;
    Scaling.unit = '';
    Scaling.Z_correction = [];
    varargout{1} = Scaling;   
    return
end

w_commands = -.68;
w_axes = .30;

this.im = [];
this.dx = [];
this.dy = [];
this.th = [];
this.xg = [];
this.yg = [];
this.tm = [];
this.i0 = [];

f = figure;
set(f,'Unit','Normalized','position',[.22 .26 .5 .5],'NumberTitle','Off',...
    'Color',[0.94 0.94 0.94])%,'IntegerHandle','Off')

uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.80+w_commands .91 .08 .04],'String','Min',...
    'horizontalalignment','center')
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.89+w_commands .91 .08 .04],'String','Max',...
    'horizontalalignment','center')

this.handles.load = uicontrol('Style','pushbutton','String', 'Load image',...
    'Units','normalized','Position',[.70+w_commands .87 .09 .08]);
this.handles.bmin = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.80+w_commands .87 .08 .05],'String','');
this.handles.bmax = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.89+w_commands .87 .08 .05],'String','');

this.handles.selectCross = uicontrol('Style','pushbutton','String', 'Select cross (dot)',...
    'Units','normalized','Position',[.70+w_commands .75 .27 .08]);
this.handles.findPoints = uicontrol('Style','slider','Value',0.8,...
    'Units','normalized','Position',[.70+w_commands .62 .27 .03]);

this.handles.unit = uicontrol('Style','edit','String', 'µm',...
    'Units','normalized','Position',[.70+w_commands .5 .08 .05]);
this.handles.pitchX = uicontrol('Style','edit','String', '1',...
    'Units','normalized','Position',[.80+w_commands .5 .08 .05]);
this.handles.pitchY = uicontrol('Style','edit','String', '1',...
    'Units','normalized','Position',[.90+w_commands .5 .08 .05]);
this.handles.fitGrid = uicontrol('Style','pushbutton','String', 'Calculate grid',...
    'Units','normalized','Position',[.70+w_commands .37 .27 .08]);
this.handles.takeVal = uicontrol('Style','pushbutton','String', 'Save',...
    'Units','normalized','Position',[.70+w_commands .02 .12 .08]);
this.handles.scalName = uicontrol('Style','edit','String', 'scaling_linear',...
    'Units','normalized','Position',[.85+w_commands .02 .12 .05]);

this.handles.dx = uicontrol('Style','edit','String', '',...
    'Units','normalized','Position',[.70+w_commands .25 .12 .05],...
    'Enable','Off');
this.handles.dy = uicontrol('Style','edit','String', '',...
    'Units','normalized','Position',[.85+w_commands .25 .12 .05],...
    'Enable','Off');
this.handles.pitchZ = uicontrol('Style','edit','String', '',...
    'Units','normalized','Position',[.70+w_commands .14 .12 .05]);
this.handles.dt = uicontrol('Style','edit','String', '',...
    'Units','normalized','Position',[.85+w_commands .14 .12 .05]);



ax = axes('position',[.04+w_axes .1 .63 .86]);
imshow(this.im,[],'Initialmagnification','fit')

this.handles.load.Callback = {@load_img,f,ax};
this.handles.selectCross.Callback = {@selectCross,f,ax};
this.handles.findPoints.Callback = {@findPoints,f,ax};
this.handles.pitchX.Callback = {@calcScaling,f};
this.handles.pitchY.Callback = {@calcScaling,f};
this.handles.fitGrid.Callback = {@fitGrid,f,ax};
this.handles.takeVal.Callback = {@takeVal,f};
this.handles.bmin.Callback = {@update_clim,f,ax};
this.handles.bmax.Callback = {@update_clim,f,ax};

setappdata(f,'this',this);


annotation(f,'textbox',[.70+w_commands .68 .27 .03],'HorizontalAlignment','center',...
    'String','Acceptance valid points','FitBoxToText','Off','EdgeColor','none');


annotation(f,'textbox',[.70+w_commands .55 .08 .05],'HorizontalAlignment','center',...
    'String','Unit','FitBoxToText','Off','EdgeColor','none');
annotation(f,'textbox',[.78+w_commands .55 .12 .05],'HorizontalAlignment','center',...
    'String','Pitch X [unit]','FitBoxToText','Off','EdgeColor','none');
annotation(f,'textbox',[.88+w_commands .55 .12 .05],'HorizontalAlignment','center',...
    'String','Pitch Y [unit]','FitBoxToText','Off','EdgeColor','none');

annotation(f,'textbox',[.68+w_commands .3 .16 .05],'HorizontalAlignment','center',...
    'String','Scaling X [unit/px]','FitBoxToText','Off','EdgeColor','none');
annotation(f,'textbox',[.83+w_commands .3 .16 .05],'HorizontalAlignment','center',...
    'String','Scaling Y [unit/px]','FitBoxToText','Off','EdgeColor','none');

annotation(f,'textbox',[.68+w_commands .19 .16 .05],'HorizontalAlignment','center',...
    'String','Z Depth [unit]','FitBoxToText','Off','EdgeColor','none');
annotation(f,'textbox',[.83+w_commands .19 .16 .05],'HorizontalAlignment','center',...
    'String','\Deltat [s]','FitBoxToText','Off','EdgeColor','none');

annotation(f,'textbox',[.83+w_commands .07 .16 .05],'HorizontalAlignment','center',...
    'String','Scaling name','FitBoxToText','Off','EdgeColor','none');

set(f,'Toolbar','figure','Menubar','None')%,'HandleVisibility','Off')
end

function load_img(~,~,f,ax)
this = getappdata(f,'this');
[filename, pathname, ok] = uigetfile({'*.tif';...
    '*.png';'*.jpg'},'Select set of images','MultiSelect','on');
if ok>0
    axes(ax), cla(ax)
    imd = imread(fullfile(pathname,filename));
    if size(imd,3)==3
        this.im = medfilt2(double(rgb2gray(imd)),[3 3]);
    else
        this.im = medfilt2(double(imd),[3 3]);
    end
    axes(ax);
    imshow(this.im,[],'Initialmagnification','fit'), hold on
    clims = get(ax,'Clim');
    set(this.handles.bmin,'string',num2str(round(clims(1))))
    set(this.handles.bmax,'string',num2str(round(clims(2))))
    setappdata(f,'this',this)
end

end

function update_clim(~,~,f,ax)
this = getappdata(f,'this');
clim1 = get(this.handles.bmin,'String');
clim2 = get(this.handles.bmax,'String');
set(ax,'Clim',[str2double(clim1) str2double(clim2)])
end

function selectCross(~,~,f,ax)

% set(f,'HandleVisibility','On');
this = getappdata(f,'this');

figure(f), axes(ax), cla(ax)
clim1 = get(this.handles.bmin,'String');
clim2 = get(this.handles.bmax,'String');
clims = [str2double(clim1) str2double(clim2)];
imshow(this.im,clims,'Initialmagnification','fit'); drawnow
htemp = imrect;
rr = round(wait(htemp)); delete(htemp)
this.tm =  this.im(rr(2):rr(2)+rr(4),rr(1):rr(1)+rr(3));
setappdata(f,'this',this)
findPoints([],[],f,ax)
% set(f,'HandleVisibility','Off')

end

function findPoints(~,~,f,ax)

this = getappdata(f,'this');

figure(f), axes(ax), cla(ax)
clim1 = get(this.handles.bmin,'String');
clim2 = get(this.handles.bmax,'String');
clims = [str2double(clim1) str2double(clim2)];
imshow(this.im,clims,'Initialmagnification','fit'); drawnow
set(f,'name','Searching points...')
drawnow

thr = get(this.handles.findPoints,'Value')/2+.5;
cc = normxcorr2(this.tm,this.im);


[xj2, yi2] = meshgrid(1:size(cc,2),1:size(cc,1));
xj2 = xj2(2:end-1,2:end-1); yi2 = yi2(2:end-1,2:end-1);

indp = cc(2:end-1,2:end-1)>cc(1:end-2,2:end-1) & cc(2:end-1,2:end-1)>cc(3:end,2:end-1) & ....
    cc(2:end-1,2:end-1)>cc(2:end-1,1:end-2) & cc(2:end-1,2:end-1)>cc(2:end-1,3:end) & ...
    cc(2:end-1,2:end-1)>thr;
xp = xj2(indp); yp = yi2(indp);
xn = xp*0; yn = yp*0;

for k = 1:length(xp)
    yn(k) = yp(k) + (log(cc(yp(k)-1,xp(k)))-log(cc(yp(k)+1,xp(k)))) /...
        ( 2*log(cc(yp(k)-1,xp(k))) - 4* log(cc(yp(k),xp(k))) + 2* log(cc(yp(k)+1,xp(k))));
    xn(k) = xp(k) + (log(cc(yp(k),xp(k)-1))-log(cc(yp(k),xp(k)+1))) /...
        ( 2*log(cc(yp(k),xp(k)-1)) - 4* log(cc(yp(k),xp(k))) + 2* log(cc(yp(k),xp(k)+1)));
end
this.xg = xn-size(this.tm,2)/2; this.yg = yn-size(this.tm,1)/2;
f1 = (xn-size(this.im,2)/2).^2 + (yn-size(this.im,1)/2).^2;
this.i0 = find(f1 == min(f1));

% set(f,'HandleVisibility','On');

plot(this.xg,this.yg,'go',this.xg,this.yg,'gx')
set(f,'name',' ')
% set(f,'HandleVisibility','Off');
setappdata(f,'this',this)
end

function fitGrid(~,~,f,ax)
this = getappdata(f,'this');
if length(this.xg)>5
    
    xf = this.xg-this.xg(this.i0); yf = this.yg-this.yg(this.i0);
    ind = abs(yf)<size(this.tm,1)/2;
    [xsort, ii] = sort(xf(ind)); ysort = yf(ind); ysort = ysort(ii);
    px = polyfit(xsort,ysort,1);
    dx0 = median(diff(xsort));
    
    ind = abs(xf)<size(this.tm,2)/2;
    [ysort, ii] = sort(yf(ind)); xsort = xf(ind); xsort = xsort(ii);
    py = polyfit(ysort,xsort,1);
    dy0 = median(diff(ysort));
    theta = atan((px(1)+py(1))/2);
    
    in0 = [dx0, dy0, theta];
    out = fminsearch(@(in) gridDist(xf,yf,in),in0);
    dx0 = out(1); dy0 = out(2); theta = out(3);
    
    limx = 0:dx0:-min(xf)+dx0; limx = [-limx(end:-1:2), 0:dx0:max(xf)+dx0];
    limy = 0:dy0:-min(yf)+dy0; limy = [-limy(end:-1:2), 0:dy0:max(yf)+dy0];
    [xdum,ydum] = meshgrid(limx,limy);
    XY = [xdum(:) ydum(:)]*[cos(theta) -sin(theta); sin(theta) cos(theta)];
    xdum(:) = XY(:,1)+this.xg(this.i0);
    ydum(:) = XY(:,2)+this.yg(this.i0);
    this.dx = dx0; this.dy = dy0; this.th = theta;
    
    %     set(f,'HandleVisibility','On');
    figure(f), axes(ax), cla(ax)
    clim1 = get(this.handles.bmin,'String');
    clim2 = get(this.handles.bmax,'String');
    clims = [str2double(clim1) str2double(clim2)];
    imshow(this.im,clims,'Initialmagnification','fit'); drawnow
    
    plot(xdum,ydum,'g',xdum',ydum','g',this.xg,this.yg,'go',this.xg,this.yg,'gx')
    %     set(f,'HandleVisibility','Off');
    setappdata(f,'this',this)
    calcScaling([],[],f)
else
    this.dx = []; this.dy = [];
    %     set(f,'HandleVisibility','On');
    figure(f), axes(ax), cla(ax)
    imshow(this.im,[],'Initialmagnification','fit'); hold on
    plot(this.xg,this.yg,'go',this.xg,this.yg,'gx')
    %     set(f,'HandleVisibility','Off');
    setappdata(f,'this',this)
end
end

function calcScaling(~,~,f)
this = getappdata(f,'this');
if ~isempty(this.dx) && ~isempty(this.dy)
    pitchX = str2double(get(this.handles.pitchX,'String'));
    pitchY = str2double(get(this.handles.pitchY,'String'));
    set(this.handles.dx,'String',num2str(pitchX*cos(this.th)/this.dx));
    set(this.handles.dy,'String',num2str(pitchY*cos(this.th)/this.dy));
else
    set(this.handles.dx,'String','');
    set(this.handles.dy,'String','');
end
setappdata(f,'this',this)
end

function takeVal(~,~,f)
this = getappdata(f,'this');

if ~isempty(this.dx) && ~isempty(this.dy) 
    pitchX = str2double(get(this.handles.pitchX,'String'));
    pitchY = str2double(get(this.handles.pitchY,'String'));
else
    disp('Missing data')
    return
end

Scaling.type = 'linear';
Scaling.X_to_x = pitchX*cos(this.th)/this.dx;
Scaling.Y_to_y = pitchY*cos(this.th)/this.dy;
Scaling.Z_to_z = eval((get(this.handles.pitchZ,'String')));
Scaling.dt = eval((get(this.handles.dt,'String')));
Scaling.unit = get(this.handles.unit,'String');
Scaling.Z_correction = [];

assignin('base',get(this.handles.scalName,'String'),Scaling)



end


function d = gridDist(xf,yf,in)
dx0 = in(1);
dy0 = in(2);
theta = in(3);

limx = 0:dx0:max(xf)+dx0; limx = [-limx(end:-1:2), limx];
limy = 0:dy0:max(yf)+dy0; limy = [-limy(end:-1:2), limy];
[xdum,ydum] = meshgrid(limx,limy);
XY = [xdum(:) ydum(:)]*[cos(theta) -sin(theta); sin(theta) cos(theta)];
xdum = XY(:,1); ydum = XY(:,2);

d = sum(min((repmat(xdum,1,length(xf)) - repmat(xf',length(xdum),1)).^2 + ...
    (repmat(ydum,1,length(yf)) - repmat(yf',length(ydum),1)).^2));
end
