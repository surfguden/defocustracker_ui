function varargout = dtracker_create(what_to_create,varargin)

switch what_to_create
    case 'preloaded_imageset' % creates an imageset of preloaded images, so you can manipulate images without writing them to HDD.
        varargout{1} = imageset_create(varargin{1});
       
    case 'dataset'
        if nargin == 10
            varargout{1} = dataset_create(varargin{1},varargin{2},varargin{3},...
                varargin{4},varargin{5},varargin{6},...
                varargin{7},varargin{8},varargin{9});
        elseif nargin==2
            varargout{1} = dataset_create(varargin{1});           
        else
            varargout{1} = dataset_create;
        end
        
    case 'imageset'
        if nargin == 3
            if ischar(varargin{2})
                fold = fullfile(varargin{1},['*.',varargin{2}]);
                D = dir(fold);
                varargout{1} = imageset_create(D(1).folder,{D(:).name});
            else
                varargout{1} = imageset_create(varargin{1},varargin{2});
            end
        elseif nargin == 2 
            fold = varargin{1};
            if fold(end-4)~='*'
                fold = fullfile(fold,'*.tif');
            end
            D = dir(fold);
            if isempty(D)
                fprintf(2, 'No image in the folder or wrong folder name. \n')
                return
            end
            varargout{1} = imageset_create(D(1).folder,{D(:).name});
        else
            varargout{1} = imageset_create;
        end
        
    case 'model'
        if nargin==1
            model_type = 'method_1';
        else
            model_type = varargin{1};
        end
        switch model_type
            case {'method_0','boundary_threshold_2d'}
                varargout{1} = process_0_create();
            case {'method_1','normalized_crosscorrelation_3d'}
                varargout{1} = calibration_1_create();
        end
        
    case 'scaling'
        if nargin==1
            scaling_type = 'linear';
        else
            scaling_type = varargin{1};
        end
        switch scaling_type
            case 'linear'
                if nargout == 0
                    scaling_linear_create()
                elseif nargout == 1
                    varargout{1} = scaling_linear_create();
                end
        end
        
    case 'tracking'
        if nargin==1
            tracking_type = 'nearest_neighbor';
        else
            tracking_type = varargin{1};
        end
        switch tracking_type
            case 'nearest_neighbor'
                varargout{1} = ptv_nearest_create();
        end
        
    otherwise
        disp(['''',what_to_create,''' not existing.'])
        return
end



