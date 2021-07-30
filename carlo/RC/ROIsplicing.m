allROI = [1:length(mask3D];
allROI = num2cell(allROI);
ROIverA = [45 42 71 52 66 67 62 32 63 56 50 31 15 61 5 8 2 9 68 11 10 7 6 80 81 13 24 23 28 79 3 17 20];
ROIverb = allROI;
for i=1:length(ROIverA)
ROIverB{1,ROIverA(1,i)} = [];
end
ROIverB = cell2mat(ROIverB);