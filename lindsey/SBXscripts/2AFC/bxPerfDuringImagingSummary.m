function area_summary = bxPerfDuringImagingSummary(mouse);
%summary of number of sessions >0.8 or 0.9 selectivity by area/layer

xlsFolder = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\';
mworksFolder = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
session_list = readtable(fullfile(xlsFolder,[mouse '_2Psessions.xlsx']));
nDays = length(session_list.date);
s = zeros(1,nDays);
b = zeros(1,nDays);
for i = 1:nDays
    filenames = dir(fullfile(mworksFolder,['data-' mouse '-' num2str(session_list.date(i)) '-*']));
    nFiles = 1:size(filenames,1);
    for ii = 1:nFiles
        load(fullfile(mworksFolder,filenames(ii).name))
        if ~input.doTestRobot & isfield(input,'trialOutcomeCell') & length(input.tGratingContrast)>20
            break
        end
    end
    trials = cell2mat(session_list.trials(i));
    if ~isempty(trials)
        str2num(trials);
    end
    [s(i), b(i)] = selectCalc(input,trials);
end

area_mat = session_list.area;
areas = unique(session_list.area);
nArea = length(areas);

depth_mat = session_list.depth;

area_summary = struct;
figure; 
[n n2] = subplotn(nArea);
for i = 1:nArea
    area_ind = find(strcmp(area_mat,areas(i)));
    subplot(n,n2,i)
    histogram(s(area_ind),[0:.1:1]); movegui('center')
    xlabel('selectivity'); ylabel('Sessions')
    ind_pt8 = find(s(area_ind)>0.8);
    ind_pt9 = find(s(area_ind)>0.9);
    title([cell2mat(areas(i)) '- N = ' num2str(length(ind_pt8)) ' > 0.8 and ' num2str(length(ind_pt9)) ' > 0.9']) 
    area_summary(i).name(1) = areas(i);
    area_summary(i).gt_pt8(1) = length(ind_pt8);
    area_summary(i).gt_pt9(1) = length(ind_pt9);
    area_summary(i).name(2) ={'L2/3'};
    area_summary(i).gt_pt8(2) = length(intersect(find(depth_mat<400),ind_pt8));
    area_summary(i).gt_pt9(2) = length(intersect(find(depth_mat<400),ind_pt9));
    area_summary(i).name(3) = {'L4'};
    area_summary(i).gt_pt8(3) = length(intersect(intersect(find(depth_mat>400),find(depth_mat<500)),ind_pt8));
    area_summary(i).gt_pt9(3) = length(intersect(intersect(find(depth_mat>400),find(depth_mat<500)),ind_pt9));
    area_summary(i).name(4) = {'L5'};
    area_summary(i).gt_pt8(4) = length(intersect(find(depth_mat>500),ind_pt8));
    area_summary(i).gt_pt9(4) = length(intersect(find(depth_mat>500),ind_pt9));
    table(area_summary(i).gt_pt8',area_summary(i).gt_pt9','RowNames',area_summary(i).name,'VariableNames',{'>0.8','>0.9'})
end



        