function d = getDistanceBetweenNeuronTypes(coordX1,coordY1,coordX2,coordY2,img,exptDate,obj,zm)
    d = cell(1,length(coordX2));
    for i = 1:length(coordX2)
        d{i} = sqrt(((coordX1-coordX2(i)).^2)+((coordY1-coordY2(i)).^2));
    end

    x_um = scalebarCalib(exptDate, obj, img, zm);

    umPerPixel = x_um./size(img,2);

    d_um = cellfun(@(x) x.*umPerPixel,d,'unif',0);
end