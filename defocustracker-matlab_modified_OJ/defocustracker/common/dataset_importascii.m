function dat = dataset_importascii(name,varargin)

if strcmp(name(end-2:end),'csv')
    T = readtable(name);
    dat = dataset_create(height(T));
    for ii = 1:length(T.Properties.VariableNames)
        myvar = table2array(T(:,ii))';
        eval(['dat.',T.Properties.VariableNames{ii},' = myvar;'])
    end
else
    
    
    A = importdata(name);
    data = A.data';
    variables = A.colheaders;
    dat = dataset_create(size(data,2));
    for ii = 1:length(variables)
        eval(['dat.',variables{ii},' = data(ii,:);'])
    end
    
    for ii = 1:size(A.textdata,1)-1
        if length(A.textdata{ii,1})>=7
            if strcmp(A.textdata{ii,1}(1:7),'scaling')
                eval(['dat.',A.textdata{ii,1},';'])
            end
        end
    end
    
end