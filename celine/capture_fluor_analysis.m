%% read in all the data
% go to WF analysis folder
close all 
clear all
cd('Z:\home\ACh\Analysis\WF');
%mice={'i2087','i2084','i2085','i2092', 'i2089'}%enter list of mice
nMice=size(mice,1);

allData=nan(6,nMice); %create an empty array to hold the data from all the mice

for iMouse = 1:nMice
    mouse=mice{iMouse}
    tableData=readtable(fullfile(['Results_', mouse,'.csv']));
    data=table2array(tableData);
    allData(1:3,iMouse)=data(:,3);
    allData(4,iMouse)=data(1,3)-data(2,3);
    allData(5,iMouse)=data(3,3)-data(2,3);
    allData(6,iMouse)=data(1,3)./data(2,3);
    clear tableData data 

end
clear mouse

