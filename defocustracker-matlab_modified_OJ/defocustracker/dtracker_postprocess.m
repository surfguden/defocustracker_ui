function varargout = dtracker_postprocess(what_to_do,varargin)
% Functions:
%   dtracker_postprocess(tracking,dataset,frame_index)
%   dtracker_postprocess('compare_true_values',dataset_meas,dataset_true,par)
%   dtracker_postprocess('merge_id',dataset_1,dataset_2)
%   dtracker_postprocess('min_cm',dataset,min_cm)
%   dtracker_postprocess('min_n_id',dataset,min_n_id)
%   dtracker_postprocess('scale',dataset)
%   dtracker_postprocess('unscale',dataset)
%   dtracker_postprocess('update_id',dataset)

% to be checked/improved:
%   dtracker_postprocess('nexscale',dataset,P)

if isstruct(what_to_do)
    tracking = what_to_do;
    what_to_do = ['tracking_',tracking.type];
end

switch what_to_do
    
    
    case 'compare_true_values'
        dat_meas = varargin{1};
        dat_true = varargin{2};
        
        if nargin==3
            par.outliers_x = 2;
            par.outliers_y = 2;
            par.outliers_z = 0.1;
            par.delta_points = 20;
        elseif nargin==4
            par = varargin{3};
        end
        if nargout==1
            errors = compare_true(dat_meas,dat_true,par);
            varargout{1} = errors;
        end
        if nargout==2
            [errors, errors_delta] = compare_true(dat_meas,dat_true,par);
            varargout{1} = errors;
            varargout{2} = errors_delta;
        end
        
%     case 'fill_empty'
%         this = varargin{1};
%         prop = {'x','y','z','dx','dy','dz','fr','id','cm'};
%         lengths = zeros(1,length(prop));
%         for ii = 1:length(prop)
%             lengths(ii) = eval(['length(this.',prop{ii},');']);
%         end
%         max_length = max(lengths);
%         for ii = 1:length(prop)
%             if lengths(ii)~=max_length
%                 eval(['this.',prop{ii},' = zeros(1,max_length);']);
%             end
%         end
%         varargout{1} = this;
        
%     case 'merge'
%         this = varargin{1};
%         pt = varargin{2};
%         varargout{1} = merge(this,pt);
%         
    case 'merge_id'
        this = varargin{1};
        pt = varargin{2};
        if ~isempty(this.id), pt.id = pt.id+max(this.id); end
        varargout{1} = merge(this,pt);
        
    case 'min_cm'
        this = varargin{1};
        min_cm = varargin{2};
        flag = this.cm>=min_cm;
        this = this(flag,:);
        varargout{1} = this;
        
    case 'min_n_id'
        this = varargin{1};
        min_id = varargin{2};
        uid = nonzeros(unique(this.id));
        flag = false(size(this.id));
        for ii = 1:length(uid)
            flag_min = this.id==uid(ii);
            if sum(flag_min)>=min_id
                flag(flag_min) = 1;
            end
        end
        this = this(flag,:);
        this = update_id(this);
        varargout{1} = this;
        
%     case 'newscale'
%         this = varargin{1};
%         P = varargin{2};
%         if this.scaling(4)==1
%             this = unscale(this);
%         end
%         this.scaling = [P.scaling_xyz, 0];
%         this = scale(this);
%         varargout{1} = this;
        
    case 'tracking_nearest_neighbor'
        this = varargin{1};
        if nargin==2
            frame_index = unique(this.fr);
        else
            frame_index = varargin{2};
        end
        this = ptv_nearest(tracking,this,frame_index);
        this.Properties.UserData.TrackingStep = tracking.tracking_step;
        varargout{1} = this;
        
    case 'scale'
        name_var = inputname(2);
        this = varargin{1};
        if nargin==3
            scal = varargin{2};
        else
            if ~isfield(this.Properties.UserData,'Scaling')
                fprintf(2, 'Scaling metadata missing. \n')
                return
            end
            scal = this.Properties.UserData.Scaling;
        end
        switch scal.type
            case 'linear'
                this = scaling_linear_apply(scal,this);
        end
        
        if nargout==1, varargout{1} = this;
        else, assignin('base',name_var,this)
        end
        
    case 'unscale'
        name_var = inputname(2);
        this = varargin{1};
        if ~isfield(this.Properties.UserData,'Scaling')
            fprintf(2, 'Scaling metadata missing \n')
            return
        end
        
        scal = this.Properties.UserData.Scaling;
        switch scal.type
            case 'linear'
                this = scaling_linear_unscale(this);
        end
        
        if nargout==1, varargout{1} = this;
        else, assignin('base',name_var,this)
        end
        
    case 'update_id'
        name_var = inputname(2);
        this = varargin{1};
        this = update_id(this);
        if nargout==1, varargout{1} = this;
        else, assignin('base',name_var,this)
        end
        
end

% function this = apply_flag(this,flag)
% if length(flag)==length(this.x)
%     prop = fields(this);
%     prop(strcmp(prop,'scaling')) = [];
%     prop(strcmp(prop,'type')) = [];
%     prop(strcmp(prop,'metadata')) = [];
%     for ii = 1:length(prop)
%         eval(['this.',prop{ii},' = this.',prop{ii},'(flag);'])
%     end
% end


% function this = scale(this)
% this.x = this.x*this.scaling(1);
% this.y = this.y*this.scaling(2);
% this.dx = this.dx*this.scaling(1);
% this.dy = this.dy*this.scaling(2);
% this.dz = this.dz*this.scaling(3);
% if this.scaling(3)>=0
%     this.z = this.z*this.scaling(3);
% else
%     this.z = -this.scaling(3)+this.z*this.scaling(3);
% end
% this.scaling(4) = 1;


% function this = unscale(this)
% this.x = this.x/this.scaling(1);
% this.y = this.y/this.scaling(2);
% this.dx = this.dx/this.scaling(1);
% this.dy = this.dy/this.scaling(2);
% this.dz = this.dz/this.scaling(3);
% if this.scaling(3)>=0
%     this.z = this.z/this.scaling(3);
% else
%     this.z = (this.scaling(3)+this.z)/this.scaling(3);
% end
% this.scaling(4) = 0;

function this = update_id(this)
iduni = nonzeros(unique(this.id));
Ntrk = length(iduni);
for kk = 1:Ntrk
    this.id(this.id==iduni(kk)) = kk;
end

% function this = merge(this,pt)
% if prod(this.scaling==pt.scaling)~=1
%     disp('Warning: Different scaling')
% end
% prop = fields(this);
% prop(strcmp(prop,'scaling')) = [];
% prop(strcmp(prop,'type')) = [];
% prop(strcmp(prop,'metadata')) = [];
% for ii = 1:length(prop)
%     eval(['this.',prop{ii},...
%         ' = [this.',prop{ii},', pt.',prop{ii},'];'])
% end
