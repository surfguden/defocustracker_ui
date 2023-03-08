function varargout = dtracker_show(var,varargin)

if istable(var)
    var_type = 'dataset';
else
    var_type = var.type;
end
name_var = inputname(1);

switch var_type
    case 'imageset'
        if nargin==1
            imageset_show(var,name_var)
        elseif nargin==2
            if isnumeric(varargin{1})
                im = imageset_read(var,varargin{1});
                varargout{1} = im;
            elseif istable(varargin{1})
                name_var2 = inputname(2);
                dataset_show(varargin{1},var,name_var2,name_var)
            end
        end
    case 'dataset'
        if nargin==1
            dataset_validate(var,name_var)
        elseif nargin>=2  
            if strcmp('plot_3d',varargin{1})
                varnames = var.Properties.VariableNames;
                varunits = var.Properties.VariableUnits;
                xyz = table2array(var(:,2:4));
                plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.'), grid on
                if strcmp(varnames{2},'X')
                    xlabel(varnames{2}), ylabel(varnames{3}), zlabel(varnames{4})
                elseif strcmp(varnames{2},'x')
                    xlabel([varnames{2},' [',varunits{2},']']) 
                    ylabel([varnames{3},' [',varunits{3},']']) 
                    zlabel([varnames{4},' [',varunits{4},']']) 
                end
                    
            elseif strcmp('plot_3d_tracks',varargin{1})
                varnames = var.Properties.VariableNames;
                varunits = var.Properties.VariableUnits;
                xyz = table2array(var(:,2:4));
                ids = unique(var.id(var.id~=0));
                if nargin==3, col = varargin{2};
                else, col = [0 0.4470 0.7410];
                end
                for ii = 1:length(ids)
                    flag = var.id==ids(ii);
                    plot3(xyz(flag,1),xyz(flag,2),xyz(flag,3),'.-',...
                        'color',col)
                    hold on
                end
                hold off
                grid on
                if strcmp(varnames{2},'X')
                    xlabel(varnames{2}), ylabel(varnames{3}), zlabel(varnames{4})
                elseif strcmp(varnames{2},'x')
                    xlabel([varnames{2},' [',varunits{2},']'])
                    ylabel([varnames{3},' [',varunits{3},']'])
                    zlabel([varnames{4},' [',varunits{4},']'])
                end
                
            elseif nargin==2 && isstruct(varargin{1})
                name_var2 = inputname(2);
                dataset_show(var,varargin{1},name_var,name_var2)
            end
        end
        
    case 'model'
        switch var.method
            case 'method_1'   
                calibration_1_show(var,name_var)
        end

end
