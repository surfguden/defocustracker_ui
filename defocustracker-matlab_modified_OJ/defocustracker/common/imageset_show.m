function varargout = imageset_show(img,varargin)

if nargin==1
    name_var = inputname(1);
elseif nargin==2
    name_var = varargin{1};
end

f = figure; set(f,'position',[50 50 800 600])

im = imageset_read(img,1,'raw');
this.n_frames = img.n_frames;
set(f,'name',['imageset: ',name_var],'numbertitle','off')

ax = axes('Parent',f,'position',[.08 .18 .9 .8]);

im2 = medfilt2(im,[3 3]);
clm = [0 (max(im2(:))+median(im2(:)))/2];

this.him = imagesc(im,clm); colormap gray
daspect([1 1 1])
axis equal
set(gca,'Color',f.Color)%,'YDir','normal')
% set(ax, 'Units', 'Normalized'); % First change to normalized units.
% set(ax, 'OuterPosition', [.08, .18, .9, .8]); % [xLeft, yBottom, width, height]
% set(gca,'PlotBoxAspectRatioMode','Auto')

this.text.stride = uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.02 .045 .1 .04],'String','Stride',...
    'horizontalalignment','center');
this.text.no_images = uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.093 .045 .28 .04],'String',['No. of images = ',num2str(this.n_frames),''],...
    'horizontalalignment','center');
this.text.min_int = uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.77 .045 .1 .04],'String','Min. intensity',...
    'horizontalalignment','center');
this.text.max_int = uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.88 .045 .1 .04],'String','Max. intensity',...
    'horizontalalignment','center');

this.img = img;

this.handles.bsk = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.02 .015 .1 .04],'String','1');
this.handles.bbw = uicontrol('Parent',f,'Style','pushbutton','units','normalized',...
    'position',[.13 .015 .05 .04],'String','<<');
this.handles.bnr = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.18 .015 .1 .04],'String','1');
this.handles.bfw = uicontrol('Parent',f,'Style','pushbutton','units','normalized',...
    'position',[.28 .015 .05 .04],'String','>>');

this.handles.bmi = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.77 .015 .1 .04],'String',num2str(clm(1)));
this.handles.bma = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.88 .015 .1 .04],'String',num2str(clm(2)));

this.handles.bnr.Callback = {@update_nr,f};
this.handles.bbw.Callback = {@update_nr_butt,f,-1};
this.handles.bfw.Callback = {@update_nr_butt,f,1};
this.handles.bmi.Callback = {@update_cl_butt,f,ax};
this.handles.bma.Callback = {@update_cl_butt,f,ax};


setappdata(f,'this',this)

if nargout==1
    varargout{1} = f;
elseif nargout==2
    varargout{1} = f;
    varargout{2} = ax;
end

end

function update_nr(~,~,f)
this = getappdata(f,'this');
n = round(str2double(this.handles.bnr.String));
if n<1, n = 1; end
if n>this.n_frames, n = this.n_frames; end
im = imageset_read(this.img,n,'raw');

set(this.him,'CData',im)
set(this.handles.bnr,'String',num2str(n))
setappdata(f,'this',this)
end

function update_nr_butt(~,~,f,ww)
this = getappdata(f,'this');
skn = round(str2double(this.handles.bsk.String))*ww;
n = round(str2double(this.handles.bnr.String))+skn;
if n<1, n = 1; end
if n>this.n_frames, n = this.n_frames; end
set(this.handles.bnr,'String',num2str(n))
setappdata(f,'this',this)
this.handles.bnr.Callback{1}([],[],f)

end

function update_cl_butt(~,~,f,ax)
this = getappdata(f,'this');
CLim1 = str2double(this.handles.bmi.String);
CLim2 = str2double(this.handles.bma.String);
if CLim2<=CLim1
    CLim2 = CLim1+1;
    set(this.handles.bma,'String',num2str(CLim2))
end
set(ax,'CLim',[CLim1 CLim2])
setappdata(f,'this',this)
end

