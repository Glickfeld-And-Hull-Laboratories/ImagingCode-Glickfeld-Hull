function [captureResults] = getCaptureValues_annulus(mice)
%getCaptureValues
%   This is a secondary function to process and combine data I extracted
%   from widefield images in ImageJ. Outputs a table with three values per
%   mouse: 1 = fluor intensity in the HTP+ region, 2 = fluor intensity in
%   the annulus, 3 = HTP+ region / annulus

nMice=size(mice,1);

allData=nan(3,nMice); %create an empty array to hold the data from all the mice

for iMouse = 1:nMice
    mouse=mice{iMouse}
    tableData=readtable(fullfile('Z:\home\ACh\Analysis\WF',['Results_', mouse,'_annulus.csv']));
    data=table2array(tableData);
    allData(1:2,iMouse)=data(:,3); %the HTP and controle means
    allData(3,iMouse)=allData(1,iMouse)/allData(2,iMouse); %HTP / control

    clear tableData data 

end
clear mouse
captureResults = allData;
end

%%change to get ratio