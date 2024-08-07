%% set path name

pn = 'G:\users\lindsey\analysisLG\081205\OriDRegistered.tif';
pn2 = 'G:\users\lindsey\dataLG\081205\081205OriD';
analpn = 'g:\users\lindsey\analysisLG\081205';

%% load stack in memory
stack_green = readtiff(pn);
stacks_green = single(stack_green);

stack_red = readtiff(pn2,[],'ch00',true);
stacks_red = single(stack_red);

%% file properties
freq = 800;
res = 512;
epoch = 64;
rep = 20;

%% plot green channel
img_green = mean(stacks_green,3);
img_green_adapt = imScale(stackLocalContrastAdaptation(img_green, 30, 1));
imagesq(img_green_adapt);

%% find cell masks
% cell find code from MH
bwimgcell = imCellEditInteractive(img_green_adapt,[]);
mask_cell = bwlabel(bwimgcell);
figure;
image(mask_cell);
figure;
imagesq(imShade(img_green_adapt,mask_cell>0));
save(fullfile(analpn,'OriDcellmasks.mat'),'mask_cell');

%% find astrocyte masks
img_red = mean(stacks_red,3);
img_red_adapt = imScale(stackLocalContrastAdaptation(img_red, 30, 1));
imagesc(img_red_adapt);
bwimgastro = imCellEditInteractive(img_red,[]);
mask_astro = bwlabel(bwimgastro);
save(fullfile(analpn,'OriDastromasks.mat'),'mask_astro')

%% reload masks
load 'g:\users\lindsey\analysisLG\081205\OriBcellmasks.mat';
load 'g:\users\lindsey\analysisLG\081205\OriBastromasks.mat';

%% remove astrocytes

common=(mask_astro>1).*mask_cell;
astro_id = unique(common(find(common(:)>0)));

mask_neuron = mask_cell;
for iCell = astro_id'
    mask_neuron(find(mask_cell==iCell))=0;
end
mask_neuron=bwlabel(mask_neuron>0);
imagesq(imShade(img_green_adapt, mask_neuron>0));
save(fullfile(analpn,'OriDneuronmasks.mat'),'mask_neuron')

%% cell time courses
load 'g:\users\lindsey\analysisLG\081205\OriDneuronmasks.mat';

%get time courses and remove low frequencies
allcells = max(unique(mask_neuron));
timeCourses = stackGetTimeCourses(stacks_green,mask_neuron);
timeCourses_lowcut = tcLowCut (timeCourses, 100, 'gaussian', 1);

%get dF/F
baseline = mean(timeCourses_lowcut(4:8:end,:));
dF = bsxfun(@minus,timeCourses_lowcut,baseline);
ratio = bsxfun(@rdivide,dF,baseline)*100;

%averaged individual cell time courses
av = tcCycleAverage(ratio,epoch);
figure;
tcOffsetPlot(tcRemoveDC(av(:,:)));
figure;
plot(av(:,:));
xlabel('Time (s)');
ylabel('Cell #');

%subplots of all averages
cycle = 1:epoch;
for iCell = 1:10;
    figure;
    plot(av(:, iCell),'r-');
    hold on
    for trial = 0:5;
        plot(ratio(cycle+(trial*epoch:trial*epoch), iCell),'k:');
        hold on;
    end
end

%% cell orientation tuning
Oris = [0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4];
Response_Means = zeros(length(Oris), allcells);
Response_Norms = zeros(length(Oris), allcells);
for iCell = 1:allcells;
    amp0 = sum(av(5:8,iCell));
    amp45 = sum(av(13:17,iCell));
    amp90 = sum(av(21:25,iCell));
    amp135 = sum(av(29:33,iCell));
    amp180 = sum(av(37:41,iCell));
    amp225 = sum(av(45:49,iCell));
    amp270 = sum(av(53:57,iCell));
    amp315 = sum(av([1 61 62 63 64],iCell));
    Response = [amp0 amp45 amp90 amp135 amp180 amp225 amp270 amp315];   
    %Response(find(Response(:) < 0)) = 0; 
    Response_Means(:,iCell) = Response;
    Response_Norms(:,iCell) = Response./max(Response);
end

%polar plots
for iCell = 11:20;
    t = Oris;
    r = Response_Norms(:,iCell)';
    figure;
    polar(t,r, '-o');
end

%circular variance for orientation
Circ_Var_Ori = zeros(1, allcells);
for n = 1:allcells; 
    Response_Mean = Response_Means(:, n);
    L = length(Oris);
    for k= 1:L
        M(k) = Response_Mean(k)*(exp(i*(2*Oris(k))));
    end
    R = (sum(M)/sum(Response_Mean));
    V = 1- abs(R);
    Circ_Var_Ori(:, n) = V;
end

%histograms
figure;
hist(Circ_Var_Ori, 20);

%CV map
mask_neuron_CV = mask_neuron;
for iCell = 1:allcells;
    mask_neuron_CV(find(mask_neuron==iCell)) = (Circ_Var_Ori(:,iCell))*100;
end
figure;
imagesc(mask_neuron_CV);
title('Circular Variance');

%CV threshold map
mask_neuron_CVthresh = mask_neuron_CV;
mask_neuron_CVthresh(find(mask_neuron_CV>91))=0;
mask_neuron_CVthresh = bwlabel(mask_neuron_CVthresh);

%Direction maps
mask_neuron_dir = mask_neuron_CVthresh;
for iCell = 1:allcells;
    mask_neuron_dir(find(mask_neuron_CVthresh==iCell)) = find(Response_Norms(:,iCell) == 1);
end
figure;
imagesc(mask_neuron_dir);
title('Direction Map for cells where CV<0.91');

%Find magnitude of response for direction opposite to preferred for measure
%of direction selectivity
DirSelectivity = [1, allcells];
for iCell = 1:allcells;
    [I J] = find(Response_Norms(:,iCell) == 1);
    if  I < 4.5;
        DirSelectivity(1,iCell) = Response_Norms(find(Response_Norms(:,iCell) == 1)+4, iCell);
    else
        DirSelectivity(1,iCell) = Response_Norms(find(Response_Norms(:,iCell)==1)-4, iCell);
    end
end

mask_neuron_dirselect = mask_neuron_CVthresh;
for iCell = 1:allcells;
    mask_neuron_dirselect(find(mask_neuron_CVthresh==iCell)) = DirSelectivity(:,iCell);
end
figure;
imagesc(mask_neuron_dirselect);
title('Direction Selectivity Map (Opp Dir/Pref Dir) for cells where CV<0.91');

mask_neuron_dirthresh = mask_neuron_dir;
mask_neuron_dirthresh(find(mask_neuron_dirselect>.6))=0;
figure;
imagesc(mask_neuron_dirthresh);
title('Direction Map for cells where CV<0.91 and DS<0.6');

%Collapse direction map into orientation map
mask_neuron_ori = mask_neuron_CVthresh;
for iCell = 1:allcells
    [I J] = find(Response_Norms(:,iCell) == 1);
    if  I > 4.5;
        mask_neuron_ori(find(mask_neuron_CVthresh==iCell)) = find(Response_Norms(:,iCell) == 1) - 4;
    else
        mask_neuron_ori(find(mask_neuron_CVthresh==iCell)) = find(Response_Norms(:,iCell) == 1);
    end
end
figure;
imagesc(mask_neuron_ori);
orimap = [1 1 1; 0 1 0; 0 0 1; 1 0 0; 1 1 0];
colormap(orimap);
title('Orientation Map for cells where CV<0.91');

