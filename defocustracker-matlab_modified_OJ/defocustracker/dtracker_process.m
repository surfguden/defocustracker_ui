function [dat, varargout] = dtracker_process(model,imageset,varargin)

if nargin==2
    frame_index = 1:imageset.n_frames;
elseif nargin==3
    frame_index = varargin{1};
end

switch model.method
    case {'method_0','boundary_threshold_2d'}
        dat = process_0(model,imageset,frame_index);
        if nargout==2
            varargout{1} = []; % process time to be implemented
        end  
        
    case {'method_1','normalized_crosscorrelation_3d'}
        [dat, process_time] = process_1(model,imageset,frame_index);
        if nargout==2
            varargout{1} = process_time;
        end    

end

if isfield(model,'tracking')
    dat = dtracker_postprocess(model.tracking,dat,frame_index);
end

if isfield(model,'scaling')
    dat = dtracker_postprocess('scale',dat,model.scaling);
end


