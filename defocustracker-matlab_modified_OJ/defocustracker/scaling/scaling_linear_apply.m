function dat = scaling_linear_apply(scaling,dat)

if ~strcmp(scaling.type,'linear')
    fprintf(2, 'Scaling not compatible \n')
    return
end

if ismember('x',dat.Properties.VariableNames)
    fprintf(2, 'Data already scaled \n')
    return
end

%%
dat.X = dat.X*scaling.X_to_x;
dat.Y = dat.Y*scaling.Y_to_y;
dat.Z = dat.Z*scaling.Z_to_z;
if scaling.Z_to_z<0
    dat.Z = -scaling.Z_to_z+dat.Z;
end
dat.DX = dat.DX*scaling.X_to_x/scaling.dt;
dat.DY = dat.DY*scaling.Y_to_y/scaling.dt;
dat.DZ = dat.DZ*scaling.Z_to_z/scaling.dt;

%%
[~, ind] = ismember('X',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'x';
dat.Properties.VariableUnits{ind} = scaling.unit;
[~, ind] = ismember('Y',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'y';
dat.Properties.VariableUnits{ind} = scaling.unit;
[~, ind] = ismember('Z',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'z';
dat.Properties.VariableUnits{ind} = scaling.unit;

[~, ind] = ismember('DX',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'vx';
dat.Properties.VariableUnits{ind} = [scaling.unit,'/s'];
[~, ind] = ismember('DY',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'vy';
dat.Properties.VariableUnits{ind} = [scaling.unit,'/s'];
[~, ind] = ismember('DZ',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'vz';
dat.Properties.VariableUnits{ind} = [scaling.unit,'/s'];

dat.Properties.UserData.Scaling = scaling;

