function model = calibration_1_create()

model.type = 'model';

model.method = 'method_1';

model.parameter.images = [];
model.parameter.images_boundary = [];
model.parameter.images_area = [];
model.parameter.z = [];
model.parameter.signal_to_noise = [];
model.parameter.similarity_map = [];
model.parameter.global_area = [];
model.parameter.n_cal = [];

model.training.imwidth = [];
model.training.imheight = [];
model.training.gauss_filter = 0;
model.training.median_filter = 5;
model.training.boundary_threshold = 3;
model.training.smoothing = 0;
model.training.n_interp_images = 0;
model.training.imresize = 1;

model.processing.scaling_xyz = [1 1 1];
model.processing.cm_final = .5;
model.processing.cm_guess = .4;
model.processing.no_overlap_radius = 25;
model.processing.subimage_interpolation = 1; % sub-image interpolation, 0 off, 1 on
model.processing.walking_step = 1; % sub-image interpolation, 0 off, 1 on
model.processing.images_in_guess_step = 9;
model.processing.n_iteration = 1;
model.processing.roi = [];



