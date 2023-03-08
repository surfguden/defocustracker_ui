function [xj, yi, cc] = find_peaks_2d(im,varargin)


if nargin==2
    thrs = varargin{1};
else
    thrs = (min(im(:))+max(im(:)))/2;
end

ims = im(2:end-1,2:end-1);
left = im(2:end-1,1:end-2)<ims;
right = im(2:end-1,3:end)<ims;
top = im(1:end-2,2:end-1)<ims;
bottom = im(3:end,2:end-1)<ims;
[II, JJ] = size(left);
kk = find(left & right & top & bottom & ims>thrs);
cc = ims(kk);
[ii, jj] = ind2sub([II,JJ],kk); yi = ii+1; xj = jj+1;

