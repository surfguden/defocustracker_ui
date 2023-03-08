function [errors, varargout] = compare_true(dat_meas,dat_true0,par)

if ~all([ismember({'fr','X','Y','Z'},dat_true0.Properties.VariableNames),...
        ismember({'fr','X','Y','Z'},dat_meas.Properties.VariableNames)])
    fprintf(2, 'Datasets are scaled or not valid for comparison \n')
    return
end

N_meas = max(dat_meas.fr);
N_true = max(dat_true0.fr);
if N_meas~=N_true
    fprintf(2, 'Warning: Frames number does not match \n')
end

disp('Error evaluation started...')

dat_true = dtracker_create('dataset',length(dat_true0.X));
dat_true.fr(:) = dat_true0.fr(:);
dat_true.X(:) = dat_true0.X(:);
dat_true.Y(:) = dat_true0.Y(:);
dat_true.Z(:) = dat_true0.Z(:);


n_measured_particles = length(dat_meas.X);
dat_meas.fr = dat_meas.fr*2-1;
dat_true.fr = dat_true.fr*2;
dat = [dat_meas; dat_true];

frame_index = 1:2:2*N_true-1;
tracking.tracking_step = 1;
tracking.bounding_box = [[-1 1]*par.outliers_x, ...
    [-1 1]*par.outliers_y, ...
    [-1 1]*par.outliers_z];

dat = ptv_nearest(tracking,dat,frame_index,'no output');

%% calculate error
flag_tp = dat.id~=0 & mod(dat.fr,2)==1;
flag_fn = dat.id==0 & mod(dat.fr,2)==0;
flag_fp = dat.id==0 & mod(dat.fr,2)==1;
tp = sum(flag_tp);
fn = sum(flag_fn);

dx = dat.DX(flag_tp);
dy = dat.DY(flag_tp);
dz = dat.DZ(flag_tp);

errors.sigma_X = sqrt(1/tp*sum(dx.^2));
errors.sigma_Y = sqrt(1/tp*sum(dy.^2));
errors.sigma_Z = sqrt(1/tp*sum(dz.^2));
errors.recall = tp/(tp+fn);
errors.precision = tp/n_measured_particles;

%% calculate error_delta

if nargout==2
    
    z_ref = dat.Z(flag_tp)+dat.DZ(flag_tp);
    z_ref_all = [z_ref; dat.Z(flag_fn)];
    z_ref_meas = [z_ref; dat.Z(flag_fp)+dat.DZ(flag_fp)];
    z = linspace(0,1,par.delta_points);
    delta_z = z(2)-z(1);
    
    errors_delta.Z = z;
    errors_delta.sigma_X = z*0;
    errors_delta.sigma_Y = z*0;
    errors_delta.sigma_Z = z*0;
    errors_delta.recall = z*0; 
    errors_delta.precision = z*0;
    errors_delta.dat = dat;
    
    for ii = 1:length(z)
        flag_tp_z = abs(z_ref-z(ii))<=delta_z;
        flag_fntp_z = abs(z_ref_all-z(ii))<=delta_z;
        flag_fptp_z = abs(z_ref_meas-z(ii))<=delta_z;
        
        errors_delta.sigma_X(ii) = sqrt(1/sum(flag_tp_z)*sum((dx(flag_tp_z)).^2));
        errors_delta.sigma_Y(ii) = sqrt(1/sum(flag_tp_z)*sum((dy(flag_tp_z)).^2));
        errors_delta.sigma_Z(ii) = sqrt(1/sum(flag_tp_z)*sum((dz(flag_tp_z)).^2));       
        errors_delta.recall(ii) = sum(flag_tp_z)/sum(flag_fntp_z);
        errors_delta.precision(ii) = sum(flag_tp_z)/sum(flag_fptp_z);
    end
    
    varargout{1} = errors_delta;
end

disp('Error evaluation done!')