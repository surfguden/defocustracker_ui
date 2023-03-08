function dat = scaling_linear_unscale(dat)


if ismember('X',dat.Properties.VariableNames) 
    fprintf(2, 'Data already unscaled. \n')
    return
end

try 
    scaling = dat.Properties.UserData.Scaling;
catch
    fprintf(2, 'Scaling not available. \n')
    return
end

dat.x = dat.x/scaling.X_to_x;
dat.y = dat.y/scaling.Y_to_y;
dat.z = dat.z/scaling.Z_to_z;
if scaling.Z_to_z<0
    dat.z = 1+dat.z/scaling.Z_to_z;
end
dat.vx = dat.vx/scaling.X_to_x*scaling.dt;
dat.vy = dat.vy/scaling.Y_to_y*scaling.dt;
dat.vz = dat.vz/scaling.Z_to_z*scaling.dt;

[~, ind] = ismember('x',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'X';
[~, ind] = ismember('y',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'Y';
[~, ind] = ismember('z',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'Z';

[~, ind] = ismember('vx',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'DX';
[~, ind] = ismember('vy',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'DY';
[~, ind] = ismember('vz',dat.Properties.VariableNames);
dat.Properties.VariableNames{ind} = 'DZ';

dat.Properties.VariableUnits = {};