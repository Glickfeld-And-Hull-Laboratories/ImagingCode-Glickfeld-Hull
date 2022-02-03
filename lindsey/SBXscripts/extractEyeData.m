function Eye_data = extractEyeData(data,rad_range)
% data is 3D movie of eye
% rad_range is [A B] where A is the min and B is the max radius
warning off;
A = cell(size(data,3),1);
B = cell(size(data,3),1);
C = cell(size(data,3),1);
D = cell(size(data,3),1);
for n = 1:size(data,3)
    A{n} = [0,0];
    B{n} = [0];
    C{n} = [0];
    D{n} = [0];
end
eye = struct('Centroid',A,'Area',B,'Val',C,'SNR',D);
radii = [];
for n = 1:size(data,3)
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.95);
              % pick the circle with best score
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
        eye(n).Val = NaN;
        eye(n).SNR = NaN;
    else
        snr = zeros(1,size(center,1));
        for idx = 1:size(center,1)
            t = double(data(:,:,n));
            vector_of_y_values = (1:size(data,1)) - center(idx,2);
            vector_of_x_values = (1:size(data,2)) - center(idx,1);
            [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
            idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
            idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
            snr(idx) = mean(t(idx1))./mean(t(idx2));
        end
        [v,idx] = max(snr);
        val = metric(idx);
        t = double(data(:,:,n));
        vector_of_y_values = (1:size(data,1)) - center(idx,2);
        vector_of_x_values = (1:size(data,2)) - center(idx,1);
        [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
        idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
        idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
        snr = mean(t(idx1))./mean(t(idx2));
        eye(n).SNR = snr;
        eye(n).Val = val;
        eye(n).Centroid = center(idx,:);
        eye(n).Area = pi*radii(idx)^2;
    end
    if mod(n,100)==0
        fprintf('Frame %d/%d\n',n,size(data,3));
    end
end
Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');
Val = double(cell2mat({eye.Val}'));
SNR = double(cell2mat({eye.SNR}'));
Eye_data.Centroid = Centroid;
Eye_data.Area = Area;
Eye_data.Val = Val;
Eye_data.SNR = SNR;

% no measurement frames
figure; 
subplot(2,2,1)
hist(sqrt(Area./pi));
xlabel('radius')
subplot(2,2,2)
hist(SNR);
xlabel('SNR')
subplot(2,2,3)
hist(Val);
xlabel('Metric')
movegui('center')

x1 = find(isnan(Area));
x2 = find(~isnan(Area));
x3 = unique([find(Val<0.1); find(Val<0.20 & SNR<1.7)]);

badFrames = unique([x1; x3]);
if length(badFrames)>25
    minx = 25;
else
    minx = length(badFrames);
end

frames = sort(randsample(length(badFrames),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,badFrames(frames(i)))); 
    hold on;
    scatter(Centroid(badFrames(frames(i)),1), Centroid(badFrames(frames(i)),2))
    title([num2str(chop(SNR(badFrames(frames(i))),2)) ' ' num2str(chop(Val(badFrames(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
movegui('center')
sgtitle(['No pupil detected- ' num2str(length(badFrames)) ' frames'])
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noPupil2.pdf']),'-dpdf','-fillpage');

goodFrames = setdiff(x2,x3);
if length(goodFrames)>25
    minx = 25;
else
    minx = length(goodFrames);
end
frames = sort(randsample(length(goodFrames),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,goodFrames(frames(i)))); 
    hold on;
    scatter(Centroid(goodFrames(frames(i)),1), Centroid(goodFrames(frames(i)),2))
    title([num2str(chop(SNR(goodFrames(frames(i))),2)) ' ' num2str(chop(Val(goodFrames(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
movegui('center')
sgtitle('Pupil detected')
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Pupil.pdf']),'-dpdf','-fillpage');

Eye_data.goodFrames = goodFrames;
Eye_data.badFrames = badFrames;
