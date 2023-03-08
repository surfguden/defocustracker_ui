function varargout = dataset_show(dat,this,varargin)

if nargin==2
    dat_name = inputname(1);
    img_name = inputname(2);
elseif nargin==4
    dat_name = varargin{1};
    img_name = varargin{2};
elseif nargin==5
    dat_name = varargin{1};
    img_name = varargin{2};
    circles = varargin{3};
end

[f, ax] = imageset_show(this);
set(f,'name',['dataset: ',dat_name,' / imageset: ',img_name],'numbertitle','off')
this = getappdata(f,'this');

if isfield(dat.Properties.UserData,'TrackingStep')
    this.tracking_step = dat.Properties.UserData.TrackingStep;
else
    this.tracking_step = 0;
end
if isfield(dat.Properties.UserData,'ParticleBoundaries')
    this.images_boundary = dat.Properties.UserData.ParticleBoundaries;
else
    this.images_boundary = [];
end

this.double_frame = 1;
if ismember('X',dat.Properties.VariableNames)
    this.dat = dat;
elseif isfield(dat.Properties.UserData,'Scaling')
    this.dat = dtracker_postprocess('unscale',dat);
else
    fprintf(2, 'Dataset not valid/Scaling missing \n')
    return
end

this.old_update = this.handles.bnr.Callback{1};
if nargin==5                        % special case for 2D-Tracking
    this.circles = circles;
end
this.handles.bnr.Callback = {@update_nr,f};



hold on
this.partc = plot(1,1,'+');
this.part = plot(1,1,'g');
this.vect = quiver(1,1,1,1,0);

if this.tracking_step>0
    btt = uicontrol('Parent',f,'Style','togglebutton',...
        'units','normalized','position',[.33 .015 .05 .04],'String','1-2');
    btt.Callback = {@next_btt,f,1};
    
    uicontrol('Parent',f,'Style','text','units','normalized',...
        'position',[.45 .045 .1 .04],'String','Vector scale',...
        'horizontalalignment','center')
    this.handles.vect_scale = uicontrol('Parent',f,'Style','edit','units','normalized',...
        'position',[.45 .015 .1 .04],'String',num2str(1));
    this.handles.vect_scale.Callback = {@update_nr,f};
    
end

if nargin==5
    % This is for the special case of the 2D-Tracking preview in order to
    % set the handle to the previously chosen frame.
    % The condition also works for the usual 2D-Tracking, as the index is
    % set to the first frame, where the circling is found.
    ind = find(~cellfun(@isempty,this.circles),1,'first');
    
    % The following adapts the starting index to the available dataset
    % offset, as the handle is chosen from the frame array, which might not
    % start with 1.
    ind = ind - (dat.fr(1)-1);
else
    ind = find(dat.id~=0,1,'first');
end

if isempty(ind), ind=1; end
set(this.handles.bnr,'String',dat.fr(ind))
setappdata(f,'this',this);
update_nr([],[],f)

if nargout==1
    varargout{1} = f;
elseif nargout==2
    varargout{1} = f;
    varargout{2} = ax;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_nr(~,~,f)
this = getappdata(f,'this');
this.old_update([],[],f)

fr = round(str2double(this.handles.bnr.String));

flag = this.dat.fr==fr;

xp = this.dat.X(flag);
yp = this.dat.Y(flag);
zp = this.dat.Z(flag);

[xf, yf] = give_positions(xp,yp,zp,this.images_boundary);

% improve this part!!
if ~isfield(this,'circles')        % don't set this one in case of 2D-Tracking
    set(this.partc,'Xdata',xp,'YData',yp)
end



set(this.part,'Xdata',xf,'YData',yf)
if ~isempty(zp) && this.tracking_step>0
    vector_scale = round(str2double(this.handles.vect_scale.String));
    dxp = this.dat.DX(flag)*vector_scale;
    dyp = this.dat.DY(flag)*vector_scale;
    set(this.vect,'Xdata',xp,'YData',yp,'UData',dxp,'VData',dyp)
else
    set(this.vect,'Xdata',0,'YData',0,'UData',0,'VData',0)
end

% try             % Case for 2D-Tracking,
%     % implements the green circles for particle boundaries.
%     if ~isempty(this.circles)
%         try
%             delete(this.circles_plot)   % resets the plot for every image
%         catch
%             disp('Try on line 133')
%         end
%         
%         % Circle Boundaries for the respective frame
%         xBound = []; yBound = [];
%         for k=1:length(this.circles{fr}.B)
%             xBound = [xBound; NaN; this.circles{fr}.B{k}(:,1)]; %#ok<AGROW>
%             yBound = [yBound; NaN; this.circles{fr}.B{k}(:,2)]; %#ok<AGROW>
%         end
%         
%         this.circles_plot = plot(xBound,yBound,"green");
%         
%         % Checks, if the underlying data is only for 2D-Tracking
%         % preview of one frame or for usual 2D-Tracking.
%         frames_encircled = size(find(~cellfun(@isempty,this.circles)),2);
%         
%         for n=1:size(this.dat.X,1)  % Here all valid particles are given
%             % their center coordinates for display
%             xp(n) = this.dat.X(n,flag);
%             yp(n) = this.dat.Y(n,flag);
%             zp(n) = this.dat.Z(n,flag);
%             this.circles_plot(n+1) = plot(xp(n),yp(n),'+');
%             str = strcat('   ID:',' ',num2str(this.dat.id(n,1)));
%             if frames_encircled > 1
%                 % In case of a preview, the IDs are not necessary.
%                 this.circles_plot(n+1+size(this.dat.X,1)) = text(xp(n),yp(n),str);
%             end
%         end
%         
%         setappdata(f,'this',this);
%     end
%     
% catch
%     disp('Try on line 130')
% end


    
end

function next_btt(~,~,f,skn)
this = getappdata(f,'this');
n = round(str2double(this.handles.bnr.String))+skn;
if n<=this.n_frames
    if this.double_frame==1
        im = imageset_read(this.img,n,'raw');
        set(this.him,'CData',im);
        this.double_frame = 2;
    else
        im = imageset_read(this.img,n-1,'raw');
        set(this.him,'CData',im);
        this.double_frame = 1;
    end
end
setappdata(f,'this',this);
end

function [xf, yf] = give_positions(x,y,z,imb)
if isempty(imb), xf = []; yf = []; return, end
z_imb = linspace(0,1,length(imb));
xf = [];
yf = [];
for ii = 1:length(x)
    [~, ind] = min(abs(z_imb-z(ii)));
    xf = [xf, NaN, x(ii)+imb{ind}(1,:)];
    yf = [yf, NaN, y(ii)+imb{ind}(2,:)];
end
end


% keyboard