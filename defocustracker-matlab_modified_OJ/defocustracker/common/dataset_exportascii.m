function dataset_exportascii(dat,name,varargin)

if nargin==3
    variables = varargin{1};
else
    variables = fields(dat);
end

% take scaling and exclude variable with not numeric values
if any(strcmp(variables,'scaling'))
    ind = find(strcmp(variables,'scaling'));
    scaling = eval(['dat.',variables{ind}]); 
    variables(ind) = [];
else
    scaling = [1 1 1 0];
end
variables(strcmp(variables,'metadata')) = [];
variables(strcmp(variables,'type')) = [];

N = zeros(1,length(variables));
for ii = 1:length(variables)
    N(ii) = eval(['length(dat.',variables{ii},')']);
end
N = max(N);

txtdat = zeros(length(variables),N);
wrong_input = false(length(variables),1);
txtformat = [];
fileID = fopen(name,'wt');

txtscaling = ['[',num2str(scaling(1),'%0.5e'),',',...
     num2str(scaling(2),'%0.5e'),',',...
     num2str(scaling(3),'%0.5e'),',',...
     num2str(scaling(4)),']'];
fprintf(fileID,'%s\n',['Defocustracker dataset, created ',date]);
fprintf(fileID,'%s%s\n','scaling = ',txtscaling);

for ii = 1:length(variables)
    var = eval(['dat.',variables{ii}]); 
    if ~isnumeric(var) || size(var,1)~=1
        wrong_input(ii) = 1;
        continue
    end   
    txtdat(ii,1:length(var)) = var;
    txtformat = [txtformat, '%12.5e '];
    fprintf(fileID,'%12s ',variables{ii});
end

fprintf(fileID,'%s\n','');

txtformat = [txtformat(1:end-1),'\n'];
txtdat(wrong_input,:) = [];
fprintf(fileID,txtformat,txtdat);
fclose(fileID);