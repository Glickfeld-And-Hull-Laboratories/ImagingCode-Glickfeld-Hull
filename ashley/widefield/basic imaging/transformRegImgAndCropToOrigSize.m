function data_reg_crop = transformRegImgAndCropToOrigSize(...
    data,regImg,movingpoints,fixedpoints)

%register data
[ypix,xpix,nfr] = size(data);
tform = fitgeotrans(movingpoints,fixedpoints,'nonreflectivesimilarity');
data_reg = imwarp(data,tform);

%crop registered image to same size as reference image
[refpoint_data,refpoint_regImg] = cpselect(data_reg(:,:,1),regImg,'Wait',true);
data_reg_crop = nan(size(data));
for i = 1:nfr
    d = data_reg(:,:,i);
    data_reg_crop(:,:,i) = d(...
        (round(refpoint_data(2))-round(refpoint_regImg(2))):...
        (round(refpoint_data(2))+(ypix-round(refpoint_regImg(2)))-1),...
        (round(refpoint_data(1))-round(refpoint_regImg(1))):...
        (round(refpoint_data(1))+(xpix-round(refpoint_regImg(1)))-1));
end
end