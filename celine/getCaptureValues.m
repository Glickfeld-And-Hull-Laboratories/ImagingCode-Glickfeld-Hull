function [captureResults] = getCaptureValues(mice)
%getCaptureValues
%   This is a secondary function to process and combine data I extracted
%   from widefield images in ImageJ. This assumes that each "results" CVS
%   has values from three ROIs: 1) the 2-photon FOV, 2) a control region,
%   and 3) the entire +HTP region, in that order. The CSVs also have to be
%   names "Results_mouseID." This script read each file by mouse ID and
%   compiles a 6 x n6Mouse array. In the array, the rows are: 1:3) raw
%   values for the three ROIs decribed above, 4) FOV - control 5) injection
%   area - control 6) FOV / control.

nMice=size(mice,1);

allData=nan(6,nMice); %create an empty array to hold the data from all the mice

for iMouse = 1:nMice
    mouse=mice{iMouse}
    tableData=readtable(fullfile('Z:\home\ACh\Analysis\WF',['Results_', mouse,'.csv']));
    data=table2array(tableData);
    allData(1:3,iMouse)=data(:,3);
    allData(4,iMouse)=data(1,3)-data(2,3);
    allData(5,iMouse)=data(3,3)-data(2,3);
    allData(6,iMouse)=data(1,3)./data(2,3);
    clear tableData data 

end
clear mouse
captureResults = allData;
end