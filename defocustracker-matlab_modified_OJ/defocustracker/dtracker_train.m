function varargout = dtracker_train(model,imset,varargin)

switch model.method
    case 'method_1'
        warning('off','all')
        if nargin==3
            model_out = calibration_1_train(model,imset,varargin{1});  
            varargout{1} = model_out;
        elseif nargin==2
            model_name = inputname(1);
            img_name = inputname(2);
            trainingset_1_create(model,imset,model_name,img_name)
        end
        warning('on','all')
end
