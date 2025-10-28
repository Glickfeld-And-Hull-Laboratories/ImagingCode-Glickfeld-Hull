mice = {'i3312'};
nMice=size(mice,1);
output_cell = cell(nMice,2);
output_cell(:,1) = mice(:,1);

% for multiple time points, naming convention should be 'Results_TIME.csv'
% TIME is the number of HOURS in NUMBER ONLY (e.g., 10 min = 0.16 hrs, so
% TIME in the file name would be 0.16)

for iMouse = 1:nMice
    allData=cell(nMice,1);
    mouse=mice{iMouse};
    fnames = ls(fullfile('G:\home\jerry\analysis\widefield\',mouse,'capture_quant','*.csv'));
    TPs = zeros(size(fnames,1),1);
    % extract time points from file names
    for iFile = 1:size(fnames,1)
        currFile = fnames(iFile,:);
        newStrs = split(currFile,'_');
        newStr = strtrim(newStrs{2}); %get rid of trailing space if exists
        currTime = newStr(1:end-4);
        TPs(iFile) = str2double(currTime);
    end
    currData = nan(4,size(TPs,1));
    extracted = nan(size(TPs,1),4);
    % find pix intensities and calculate ratio
    for iTP = 1:size(TPs)
        time = TPs(iTP);
        tableData=readtable(fullfile('G:\home\jerry\analysis\widefield\',mouse,'capture_quant',['Results_' num2str(time) '.csv']));
        data=table2array(tableData);
        currData(1:2,iTP)=data(:,3); %the HTP and controle means
        currData(3,iTP)=currData(1,iTP)/currData(2,iTP); %HTP / control
        currData(4,iTP) = time;
    end
    extracted(:,1) = currData(3,:);
    extracted(:,2) = currData(4,:);
    extracted(:,3) = currData(1,:); % HTP
    extracted(:,4) = currData(2,:); % control
    sorted_vals = sortrows(extracted,2);
    output_cell{iMouse,2} = sorted_vals;
    %clear tableData data 
end

captureResults = output_cell;

warning('on');