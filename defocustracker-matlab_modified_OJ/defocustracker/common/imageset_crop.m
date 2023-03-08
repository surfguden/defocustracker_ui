function IMc = imageset_crop(IM0,X,Y,imrect)

% imrect = [width height];
IW = flip(imrect);
xj = (-(IW(2)-1)/2:(IW(2)-1)/2)+X;
yi = (-(IW(1)-1)/2:(IW(1)-1)/2)+Y;

dd = zeros(1,4); xx = dd; yy = dd;
IMs = cell(1,4);


if floor(xj(1))<1 || floor(yi(1))<1 ||...
        ceil(xj(end))>size(IM0,2) || ceil(yi(end))>size(IM0,1)
    IMc = []; return
end

% if mod(X,1)==0 && mod(Y,1)==0
%     keyboard
%     IMc = IM0(yi,xj); return
% end

IMs{1} = IM0(floor(yi),floor(xj));
xx(1) = xj(1)-floor(xj(1)); yy(1) = yi(1)-floor(yi(1));
IMs{2} = IM0(ceil(yi),floor(xj));
xx(2) = xj(1)-floor(xj(1)); yy(2) = yi(1)-ceil(yi(1));
IMs{3} = IM0(floor(yi),ceil(xj));
xx(3) = xj(1)-ceil(xj(1)); yy(3) = yi(1)-floor(yi(1));
IMs{4} = IM0(ceil(yi),ceil(xj));
xx(4) = xj(1)-ceil(xj(1)); yy(4) = yi(1)-ceil(yi(1));
dd = sqrt(xx.^2+yy.^2);
[~, ii] = sort(dd);
xx = abs(xx(ii(1:3))); yy = abs(yy(ii(1:3)));
a = -xx(1)-yy(1)+1; b = xx(1); c = yy(1);
IMc = IMs{ii(1)}*a+IMs{ii(2)}*b+IMs{ii(3)}*c;

