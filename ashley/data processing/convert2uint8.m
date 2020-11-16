function img8 = convert2uint8(img)

img8 = uint8((img./max(img(:))).*255);

end