function trainingset_1_create(model,img,varargin)


if nargin==4
    model_name = varargin{1};
    img_name = varargin{2};
elseif nargin==2
    model_name = inputname(1);
    img_name = inputname(2);
end
    
stretch_factor = 5/4;

dat = dtracker_create('dataset',img.n_frames);
[f, ax] = dataset_show(dat,img);
f.Position(3) = f.Position(3)*stretch_factor;
ax.Position([1 3]) = ax.Position([1 3])/stretch_factor;
set(f,'name',['training set on: ',img_name,';  model: ',model_name],'numbertitle','off')

this = getappdata(f,'this');
this.model = model;
hfields = fields(this.handles);
for ii = 1:length(hfields)
    eval(['this.handles.',hfields{ii},'.Position([1 3]) = ',...
        'this.handles.',hfields{ii},'.Position([1 3])/stretch_factor;'])
end
tfields = fields(this.text);
for ii = 1:length(tfields)
    eval(['this.text.',tfields{ii},'.Position([1 3]) = ',...
        'this.text.',tfields{ii},'.Position([1 3])/stretch_factor;'])   
end

edit_width = this.handles.bmi.Position(3);

block_height = .93;
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.8 block_height .18 .04],'String','Particle position (x,y)',...
    'horizontalalignment','center')
this.handles.xpos = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.8 block_height-.03 edit_width .04],'String',num2str(0));
this.handles.ypos = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.98-edit_width block_height-.03 edit_width .04],'String',num2str(0));

block_height = .84;
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.8 block_height .18 .04],'String','Bounding box size (width, height)',...
    'horizontalalignment','center')
this.handles.imwidth = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.8 block_height-.03 edit_width .04],'String',num2str(0));
this.handles.imheight = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.98-edit_width block_height-.03 edit_width .04],'String',num2str(0));

block_height = .74;
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.8 block_height .09 .04],'String','First frame',...
    'horizontalalignment','center')
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.89 block_height .09 .04],'String','Last frame',...
    'horizontalalignment','center')
this.handles.fr_st = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.8 block_height-.03 edit_width .04],'String',num2str(1));
this.handles.fr_end = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.98-edit_width block_height-.03 edit_width .04],'String',num2str(img.n_frames));
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.8 block_height-.08 .09 .04],'String','Stride',...
    'horizontalalignment','center')
this.handles.fr_stride = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.8 block_height-.11 edit_width .04],'String',num2str(1));



this.handles.apply = uicontrol('Parent',f,'Style','pushbutton','units','normalized',...
    'position',[.8 .52 .18 .08],'String','Apply');


block_height = .4;
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.8 block_height .09 .04],'String','Median filter',...
    'horizontalalignment','center')
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.88 block_height .11 .04],'String','Bounding threshold',...
    'horizontalalignment','center')
this.handles.median_filter = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.8 block_height-.03 edit_width .04],'String',...
    num2str(model.training.median_filter));
this.handles.boundary_threshold = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.98-edit_width block_height-.03 edit_width .04],'String',...
    num2str(model.training.boundary_threshold));
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.8 block_height-.08 .09 .04],'String','Smoothing',...
    'horizontalalignment','center')
uicontrol('Parent',f,'Style','text','units','normalized',...
    'position',[.89 block_height-.08 .09 .04],'String','No. interp images',...
    'horizontalalignment','center')
this.handles.smoothing = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.8 block_height-.11 edit_width .04],'String',...
    num2str(model.training.smoothing));
this.handles.n_interp_images = uicontrol('Parent',f,'Style','edit','units','normalized',...
    'position',[.98-edit_width block_height-.11 edit_width .04],'String',...
    num2str(model.training.n_interp_images));







this.handles.train = uicontrol('Parent',f,'Style','pushbutton','units','normalized',...
    'position',[.8 .18 .18 .08],'String','Train');


this.handles.apply.Callback = {@apply,f};
this.handles.train.Callback = {@train,f,model_name};

setappdata(f,'this',this)


function apply(~,~,f)
this = getappdata(f,'this');
fr_st = str2double(get(this.handles.fr_st,'String'));
fr_end = str2double(get(this.handles.fr_end,'String'));
fr_stride = str2double(get(this.handles.fr_stride,'String'));
ind = fr_st:fr_stride:fr_end;
this.dat = dtracker_create('dataset',length(ind));
this.dat.fr = ind;
this.dat.x(:) = str2double(get(this.handles.xpos,'String'));
this.dat.y(:) = str2double(get(this.handles.ypos,'String'));
this.dat.z = (ind-ind(1))/(ind(end)-ind(1));
imw = str2double(get(this.handles.imwidth,'String'));
imh = str2double(get(this.handles.imheight,'String'));
this.images_boundary{1} = [-imw -imw imw imw -imw;...
     -imh imh imh -imh -imh]/2;

setappdata(f,'this',this)
this.handles.bnr.Callback{1}([],[],f)

function train(~,~,f,model_name)
this = getappdata(f,'this');
%%
model = this.model;
model.training.imwidth = str2double(get(this.handles.imwidth,'String'));
model.training.imheight = str2double(get(this.handles.imheight,'String'));
model.training.median_filter = str2double(get(...
    this.handles.median_filter,'String'));
model.training.boundary_threshold = str2double(get(...
    this.handles.boundary_threshold,'String'));
model.training.smoothing = str2double(get(...
    this.handles.smoothing,'String'));
model.training.n_interp_images = str2double(get(...
    this.handles.n_interp_images,'String'));

dat = this.dat;
img = this.img;

set(this.handles.train,'String','Training...')
drawnow
model = dtracker_train(model,img,dat);
set(this.handles.train,'String','Train')

eval([model_name,' = model;'])
eval(['dtracker_show(',model_name,')'])


assignin('base',model_name,model)
setappdata(f,'this',this)


