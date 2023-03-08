function [dat, process_time] = process_1(model,imageset,frame_index)

disp('Evaluation started ...')

cal = model.parameter;
cal.imwidth = model.training.imwidth;
cal.gauss_imheight = model.training.imheight;
cal.gauss_filter = model.training.gauss_filter;
cal.median_filter = model.training.median_filter;
cal.imresize = model.training.imresize;

par = model.processing;

dat = dtracker_create('dataset');
eval_times = frame_index*0;

ind = find(frame_index<=imageset.n_frames,1,'last');
frame_index = frame_index(1:ind);

n_frames = length(frame_index);

parfor n = 1:n_frames
    im = imageset_read(imageset,(frame_index(n)));
    
    if ~isempty(par.roi)
        im = imcrop(im, par.roi);
    end
  
     
    [datd, eval_times(n)] = singleprocess_1(im,cal,par);
    datd.fr(:) = frame_index(n);
    dat = [dat; datd];
    timestr = secs2hms(eval_times(n)*(n_frames-n));
    if n==n_frames, timestr = ' '; end
    
    note_eval = [num2str(n),'/',num2str(n_frames),' - ',...
        num2str(round(eval_times(n))),' s - ',...
        num2str(length(datd.X)),' particles - time left: ',...
        timestr];
    
    disp(note_eval)
end

if ~isempty(par.roi)
    dat.X = dat.X+par.roi(1);
    dat.Y = dat.Y+par.roi(2);
end

dat.Properties.UserData.ParticleBoundaries = model.parameter.images_boundary;


process_time = sum(eval_times);
timestr = secs2hms(process_time);
disp('----------------------------')
disp(['Evaluation done! Total Time: ',timestr])




