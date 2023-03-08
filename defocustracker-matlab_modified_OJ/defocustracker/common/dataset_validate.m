function dataset_validate(dat,varargin)

if nargin==2
    name_var = varargin{1};
else
    name_var = inputname(1);
end


f = figure;
set(f,'name',['dataset: ',name_var],'numbertitle','off')
ax = axes('position',[.06 .14 .64 .8]);

x = table2array(dat(:,2));
y = table2array(dat(:,3));

hpt.blue = plot(x,y,'b.'); grid on, hold on
hpt.green = plot(x(1),y(1),'g.');
hpt.grey = plot(x(1),(1),'.','color',[.6 .6 .6]);
xlim([min(x) max(x)]), ylim([min(y) max(y)])

flag = true(size(x));
flag_undo = [];
flag_tag = false(size(x));

setappdata(f,'flag',flag);
setappdata(f,'name_var',name_var);
setappdata(f,'flag_undo',flag_undo);
setappdata(f,'flag_tag',flag_tag);
setappdata(f,'dat',dat);

vars = fieldnames(dat); vars = vars(1:9);

% var1 drop down
annotation('textbox',[.06 .04 .03 .04],...
    'String','X: ','FitBoxToText','off','EdgeColor','none');
var1 = uicontrol('Style', 'popup','string',vars,...
    'Units','normalized','Position',[.09 -.02 .16 .1]);
set(var1,'Value',2)

% var2 drop down
annotation('textbox',[.3 .04 .03 .04],...
    'String','Y: ','FitBoxToText','off','EdgeColor','none');
var2 = uicontrol('Style', 'popup','string',vars,...
    'Units','normalized','Position',[.33 -.02 .16 .1]);
set(var2,'Value',3)

% buttons remove inside
rin = uicontrol('Style', 'pushbutton','String', 'Remove Inside',...
    'Units','normalized','Position',[.73 .72 .12 .1]);

% buttons remove outside
rout = uicontrol('Style', 'pushbutton','String', 'Remove Outside',...
    'Units','normalized','Position',[.73 .84 .12 .1]);

% button tag
tag = uicontrol('Style', 'pushbutton','String', 'Tag',...
    'Units','normalized','Position',[.86 .84 .12 .1]);

% button clear tag
notag = uicontrol('Style', 'pushbutton','String', 'Clear Tags',...
    'Units','normalized','Position',[.86 .72 .12 .1]);


% button undo
und = uicontrol('Style', 'pushbutton','String', 'Undo',...
    'Units','normalized','Position',[.73 .55 .12 .09]);

% button reset
res = uicontrol('Style', 'pushbutton','String', 'Reset',...
    'Units','normalized','Position',[.86 .55 .12 .09]);

% edit marker size
annotation('textbox',[.73 .29 .13 .04],...
    'String','Marker size','FitBoxToText','Off','EdgeColor','none');
mrk = uicontrol('Style', 'edit','String','6','Units','normalized',...
    'Position',[.73 .25 .13 .04]);

% check boxes rectangle/polygonal
bg = uibuttongroup('Units','normalized','Position',[.73 .4 .25 .05]);
rec = uicontrol(bg,'Style','radiobutton','String','Rectangle',...
    'Units','normalized','Position',[0 0 .5 1]);
uicontrol(bg,'Style','radiobutton','String','Polygonal',...
    'Units','normalized','Position',[0.5 0 .5 1]);

% apply
app = uicontrol('Style', 'pushbutton','String', 'Save',...
    'Units','normalized','Position',[.73 .04 .12 .1]);

% save as
sas = uicontrol('Style', 'pushbutton','String', 'Save as',...
    'Units','normalized','Position',[.86 .04 .12 .1]);


var1.Callback = {@update_var,f,hpt,ax,var1,var2};
var2.Callback = {@update_var,f,hpt,ax,var1,var2};
mrk.Callback = {@update_mrk,hpt};
rin.Callback = {@remove,f,hpt,ax,rec,1,var1,var2};
rout.Callback = {@remove,f,hpt,ax,rec,0,var1,var2};
tag.Callback = {@make_tag,f,hpt,ax,rec,var1,var2};
notag.Callback = {@clear_tag,f,hpt,var1,var2};
res.Callback = {@reset,f,hpt,var1,var2};
und.Callback = {@undo,f,hpt,var1,var2};
app.Callback = {@apply,f};
sas.Callback = {@save_as,f};


title(['Total particles: ',num2str(length(x))])

set(f,'position',[200 100 800 600],'Toolbar','figure','Menubar','None','NumberTitle','Off')
update(f,hpt,var1, var2)


function update(f,hpt,var1,var2)
flag = getappdata(f,'flag');
flag_tag = getappdata(f,'flag_tag');
[x,y] = get_x_y(f,var1,var2);
set(hpt.blue,'XData',x(flag),'YData',y(flag))
set(hpt.grey,'XData',x(~flag),'YData',y(~flag))
set(hpt.green,'XData',x(flag_tag & flag),'YData',y(flag_tag & flag))


function update_var(~,~,f,hpt,ax,var1,var2)
update(f,hpt,var1,var2)
[x,y] = get_x_y(f,var1,var2);
set(ax,'Xlim',[min(x) max(x)+.00001],'YLim',[min(y) max(y)+.00001])


function update_mrk(~,event,hpt)
mrk_size = str2double(event.Source.String);
set(hpt.blue,'MarkerSize',mrk_size)
set(hpt.grey,'MarkerSize',mrk_size)
set(hpt.green,'MarkerSize',mrk_size)


function remove(~,~,f,hpt,ax,rec,is_in,var1,var2)

flag_tag = getappdata(f,'flag_tag');
flag_old = getappdata(f,'flag');
[x,y] = get_x_y(f,var1,var2);

switch rec.Value
    case 0
        hrc = impoly(ax); vxy = wait(hrc); delete(hrc)
        flag = inpolygon(x,y,vxy(:,1),vxy(:,2));
    case 1
        hrc = imrect(ax); vxy = wait(hrc); delete(hrc)
        flag = x>vxy(1) & y>vxy(2) & x<vxy(1)+vxy(3) & y<vxy(2)+vxy(4);
end

if sum(flag_tag)~=0
    if is_in==1, flag = ~(flag & flag_tag);
    else, flag = ~(~flag & flag_tag);
    end
elseif is_in==1
    flag = ~flag;
end

flag(~flag_old) = 0;
setappdata(f,'flag_undo',flag_old)
setappdata(f,'flag',flag)
update(f,hpt,var1,var2)


function make_tag(~,~,f,hpt,ax,rec,var1,var2)

[x,y] = get_x_y(f,var1,var2);

switch rec.Value
    case 0
        hrc = impoly(ax); vxy = wait(hrc); delete(hrc)
        flag = inpolygon(x,y,vxy(:,1),vxy(:,2));
    case 1
        hrc = imrect(ax); vxy = wait(hrc); delete(hrc)
        flag = x>vxy(1) & y>vxy(2) & x<vxy(1)+vxy(3) & y<vxy(2)+vxy(4);
end
setappdata(f,'flag_tag',flag)
update(f,hpt,var1,var2)

function reset(~,~,f,hpt,var1,var2)
flag = getappdata(f,'flag');
flag(:) = true;
setappdata(f,'flag',flag)
setappdata(f,'flag_undo',flag)
update(f,hpt,var1,var2)

function undo(~,~,f,hpt,var1,var2)
flag_undo = getappdata(f,'flag_undo');
if ~isempty(flag_undo)
    setappdata(f,'flag',flag_undo)
    setappdata(f,'flag_undo',[])
    update(f,hpt,var1,var2)
end

function clear_tag(~,~,f,hpt,var1,var2)
flag = getappdata(f,'flag_tag');
flag(:) = false;
setappdata(f,'flag_tag',flag)
update(f,hpt,var1,var2)

function apply(~,~,f)

flag = getappdata(f,'flag');
dat = getappdata(f,'dat');
name_var = getappdata(f,'name_var');
dat = dtracker_postprocess('apply_flag',dat,flag);
assignin('base',name_var,dat)
close(f)

function save_as(~,~,f)
name_var = getappdata(f,'name_var');
titleInputdlg = 'New data name';
nameIsValid = 0; %existAlready = 1;
while nameIsValid==0 %~(existAlready==0 && nameIsValid==1)
    answer = inputdlg(titleInputdlg,'',1,{name_var});
    if isempty(answer), return, end
    % existAlready = exist(answer{1},'var');
    nameIsValid = isvarname(answer{1});
    titleInputdlg = 'Name not valid, please change';
end

flag = getappdata(f,'flag');
dat = getappdata(f,'dat');
dat = dat(flag,:);
assignin('base',answer{1},dat)
close(f)


function [x,y] = get_x_y(f,var1,var2)
dat = getappdata(f,'dat');
v1 = var1.Value;
v2 = var2.Value;
x = eval(['dat.',var1.String{v1},';']);
y = eval(['dat.',var2.String{v2},';']);