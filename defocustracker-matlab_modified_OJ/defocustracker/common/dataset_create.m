function this = dataset_create(varargin)

var = {'fr','X','Y','Z','DX','DY','DZ','id','cm'};


if nargin==4 || nargin==8 || nargin==9
    A = [];
    for ii = 1:nargin
        A = [A, varargin{ii}(:)];
    end
    this = array2table(A);
    this.Properties.VariableNames = var(1:nargin);  
elseif nargin==1
    A = zeros(varargin{1},9);
    this = array2table(A);
    this.Properties.VariableNames = var;      
else
    A = zeros(1,9);
    this = array2table(A);
    this.Properties.VariableNames = var;
    this(:,:) = [];
end
% 
% this.scaling = [1 1 1 0];
% this.type = 'dataset';
% this.metadata.images_boundary = [];
% this.metadata.tracking_step = 0;

