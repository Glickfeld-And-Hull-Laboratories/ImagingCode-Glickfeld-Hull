a = load('Z:\home\tj\2p_Imaging\i2566\230908_i2566_expression\snapshot_2.0zoom_15power_250depth.mat')
a = a.img; 
green_bright_a = imlocalbrighten(a, 0.1)
figure; imagesc(green_bright_a)

a_gray = rgb2gray(a)
figure; imshow(a_gray)
a_bright = imlocalbrighten(a_gray, 0.2)

figure;
montage({a_gray,a_bright})


%%
raw_green = load('Z:\home\tj\2p_Imaging\i2567\230908_i2567_expression\snapshot_8.0zoom_20power_231depth.mat');
raw_green = raw_green.img;
figure; imagesc(raw_green);

green_bright = imlocalbrighten(raw_green,0.2);
figure; imagesc(green_bright);

raw_red = load('Z:\home\tj\2p_Imaging\i2567\230908_i2567_expression\snapshot_red_8.0zoom_40power_231depth.mat');
raw_red = raw_red.img;
figure; imagesc(raw_red);

red_bright = imlocalbrighten(raw_red,0.2);
figure; imagesc(red_bright);

green_gray = rgb2gray(raw_green);
figure; imshow(green_gray);

green_gray_bright = imlocalbrighten(green_gray, 0.8);
figure; imshow(green_gray_bright);
