function C = normxcorr2_general(im,imt,iswin)

if iswin
    C = normxcorr2_mex(im,imt);
else
    C = normxcorr2(im,imt);
end