function mydat = ptv_nearest(tracking,mydat,frame_index,varargin)
%
% tracking.tracking_step is the frame distance to the second frame
% tracking.bounding_box is [minDX maxDX minDY maxDY minDZ maxDZ]
% dat is a DefocusTracker dataset
% frame_index is a vector with frames index
%
% The ptv analysis is performed with the following scheme:
% (frame_index(1)+tracking_step) - (frame_index(1))
% (frame_index(2)+tracking_step) - (frame_index(2))
% (frame_index(3)+tracking_step) - (frame_index(3))
% .....


% dat = dtracker_postprocess('unscale',dat);

fr_1 = frame_index;
fr_2 = fr_1+tracking.tracking_step;
indgood = fr_2<=max(mydat.fr);
fr_1 = fr_1(indgood);
fr_2 = fr_2(indgood);

par.minmax_DX = tracking.bounding_box(1:2);
par.minmax_DY = tracking.bounding_box(3:4);
par.minmax_DZ = tracking.bounding_box(5:6);


if length(fr_1)==1
    fr_step = tracking.tracking_step;
else
    fr_step = fr_1(2)-fr_1(1);
end

if nargin==3
    disp('Tracking - nearest neighbor - started...')
end
tic

dat.fr = mydat.fr';
xyz = table2array(mydat(:,2:4));
dat.x = xyz(:,1)';
dat.y = xyz(:,2)';
dat.z = xyz(:,3)';
dat.dx = dat.x*0;
dat.dy = dat.x*0;
dat.dz = dat.x*0;
dat.id = dat.x*0;

if fr_step==tracking.tracking_step
    lastID = 0;
    for n = 1:length(fr_1)
        [dat, id1, id2]  = ptv(par,dat,fr_1(n),fr_2(n));
        id10 = find(dat.id(id1)==0);
        newID = lastID+(1:length(id10));
        dat.id(id1(id10)) = newID;
        dat.id(id2) = dat.id(id1);
        lastID = max(dat.id);
    end
else
    for n = 1:length(fr_1)
        [dat, id1, id2] = ptv(par,dat,fr_1(n),fr_2(n));
        dat.id(id1) = 1;
        dat.id(id2) = 2;
    end
end

mydat(:,2:8) = array2table([dat.x; dat.y; dat.z; dat.dx; dat.dy; dat.dz; dat.id]');

mydat.Properties.UserData.TrainingStep = tracking.tracking_step;
totTime = toc;
minTime = floor(totTime/60);
secTime = round(totTime-minTime*60);


if nargin==3
    disp(['Tracking done! Total time: ',num2str(minTime),' min ',num2str(secTime),' sec'])
end

function [ell1, varargout] = ptv(par,ell1,n_1,n_2)

ind01 = find(ell1.fr==n_1);
ind02 = find(ell1.fr==n_2);

if ~isempty(ind01) && ~isempty(ind02)
    
    x1 = ell1.x(ind01); y1 = ell1.y(ind01); z1 = ell1.z(ind01);
    x2 = ell1.x(ind02); y2 = ell1.y(ind02); z2 = ell1.z(ind02);
    %     ell1.dx(ind01) = 0; ell1.dy(ind01) = 0; ell1.dz(ind01) = 0;
    %     ell1.dx(ind02) = 0; ell1.dy(ind02) = 0; ell1.dz(ind02) = 0;
    ind1 = []; ind2 = [];
    
    d1 = length(x1);
    d2 = length(x2);
    ind_d1 = 1:d1;
    ind_d2 = 1:d2;
    
    diff_x = repmat(x2, [d1 1]) - repmat(x1', [1 d2]);
    diff_y = repmat(y2, [d1 1]) - repmat(y1', [1 d2]);
    diff_z = repmat(z2, [d1 1]) - repmat(z1', [1 d2]);
    
    mask = false(size(diff_x));
    mask(diff_x<par.minmax_DX(1) | diff_x>par.minmax_DX(2))=1;
    mask(diff_y<par.minmax_DY(1) | diff_y>par.minmax_DY(2))=1;
    mask(diff_z<par.minmax_DZ(1) | diff_z>par.minmax_DZ(2))=1;
    
    diff_xyz = sqrt(diff_x.^2+diff_y.^2);%+diff_z.^2);
    diff_xyz(mask) = NaN;
    
    [~, jj] = find(diff_xyz==min(diff_xyz(:)),1);
    while ~isempty(jj) && ~isempty(diff_xyz)
        [ii, jj] = find(diff_xyz==min(diff_xyz(:)),1);
        ind1 = [ind1 ind_d1(ii)];
        ind2 = [ind2 ind_d2(jj)];
        diff_xyz(ii,:) = []; diff_xyz(:,jj) = [];
        ind_d1(ii) = []; ind_d2(jj) = [];
    end
    
    if ~isempty(ind2)
        ell1.dx(ind01(ind1)) = x2(ind2)-x1(ind1);
        ell1.dy(ind01(ind1)) = y2(ind2)-y1(ind1);
        ell1.dz(ind01(ind1)) = z2(ind2)-z1(ind1);
        ell1.dx(ind02(ind2)) = x2(ind2)-x1(ind1);
        ell1.dy(ind02(ind2)) = y2(ind2)-y1(ind1);
        ell1.dz(ind02(ind2)) = z2(ind2)-z1(ind1);
    end
    
    if nargout == 3
        varargout{1} = ind01(ind1);
        varargout{2} = ind02(ind2);
    end
    
elseif nargout == 3
    varargout{1} = [];
    varargout{2} = [];
end

