clearvars
close all


%%
dataset = 'oriAdapt_V1';
eval(dataset);

i475_expts = 30:79;
i472_expts = 80:125;


exp_list = 30:125;

data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';

%% 

for i = 1:length(expt)
  %% Import Session info 
%     iexp = exp_list(i);
    
    mouse = expt(i).mouse;
    date = expt(i).date;

    nrun = size(expt(i).runs,1);
    
    run_str = ['runs']; 
    run_str = [run_str '-' expt(i).runs(1,:)];

    % load initial mask info to be used later in loop
%     load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

 %%  Get Input Struct

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
    if strcmp(expt(i).folder,'lindsey')
        data_base = LG_base;
    elseif strcmp(expt(i).folder,'camaron')
        data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
    end
 
 
    CD = [data_base '\Data\2P_images\' expt(i).mouse '\' expt(i).date '\' expt(i).runs(1,:)];
    cd(CD);
    imgMatFile = [expt(i).runs(1,:) expt(i).runs_suffix(1,:) '.mat']; % DONE; Make variable to pull for oriAdapt_V1 that points to imgMatFile of restarted runs (ex: 001_000_001)
    load(imgMatFile);
%     fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).time_mat(1,:) '.mat'];
%     load(fName);


    
    z_pos(i) = info.config.knobby.pos.z; 
end
