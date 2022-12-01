%% Get good TC extractions that meet conditions such as: behaving, passive, direction tuning, 
% FB OFF, and s > 0.8 (can modify to pull experiments that meet other conditions of interest 
% based on oriAdapt_V1_cam and input.mat)

%%
dataset = 'oriAdapt_V1_cam';
eval(dataset);

i1402_expts = 1:16;
i1403_expts = 17:29;
i475_expts = 30:79;
i472_expts = 80:131;
i1369_expts = 137;
i1370_expts = 138;
i1372_expts = 139;
i484_expt = []; % TC's need to be extracted

% n_i475_all = len(i475_expts);
% n_i472_all = len(i472_expts);
    
%% Days with passive run and direction tuning
% i475_pd_run_list = [];
% i472_pd_run_list = [];
% 
% for i = i475_expts
%     if ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         i475_pd_run_list = [i475_pd_run_list i];
%     end
% end
% 
% for i = i472_expts
%     if ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         i472_pd_run_list = [i472_pd_run_list i];
%     end
% end
% 
% 
% n_i475_pd_runs = len(i475_pd_run_list)
% n_i472_pd_runs = len(i472_pd_run_list)


%% Get good extractions with passive and direction tuning

% i475_good_list = [];
% i472_good_list = [];
% 
% for i = i475_expts
%     if expt(i).TCs_extracted == 1 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         i475_good_list = [i475_good_list i];
%     end
% end
% 
% for i = i472_expts
%     if expt(i).TCs_extracted == 1 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         i472_good_list = [i472_good_list i];
%     end
% end
% 
% % n_i475_done = len(i475_good_list);
% % n_i472_done = len(i472_good_list);
% 
% % i475_prcnt_comp = n_i475_done/n_i475_pd_runs
% % i472_prcnt_comp = n_i472_done/n_i472_pd_runs
% 
% % i475_good_list(1:5) = [];
% % expts_to_update = [i475_good_list i472_good_list];


%% Get incomplete days with passive and direction tuning

% i475_incomplete = [];
% i472_incomplete = [];
% 
% 
% for i = i475_expts
%     if expt(i).TCs_extracted == 0 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         i475_incomplete = [i475_incomplete i];
%     end
% end
% 
% for i = i472_expts
%     if expt(i).TCs_extracted == 0 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 
%         i472_incomplete = [i472_incomplete i];
%     end
% end
% 
% disp(i475_incomplete)
% disp(i472_incomplete)


%%

i1402_good_list = genExptListNoPassive(i1402_expts);
%%

i1403_good_list = genExptListNoPassive(i1403_expts);
%%

i475_good_list = [];

for i = i475_expts
    if expt(i).TCs_extracted == 1 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 

        irun = 1;
        nrun = 1;
        run_str = ['runs']; 
        run_str = [run_str '-' expt(i).runs(irun,:)];

        LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        if strcmp(expt(i).folder,'lindsey')
            data_base = LG_base;
        elseif strcmp(expt(i).folder,'camaron')
            data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
        end

        mouse = expt(i).mouse;
        date = expt(i).date;
    
        load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
        if any(celleqel2mat_padded(adapt_input.tDoFeedbackMotion(21:end)))
            continue
        else
            ntrials = length(adapt_input.tDoFeedbackMotion);
            [s b] = selectCalc(adapt_input,21:ntrials);
            if s>=0.8 & abs(b)<0.1
                i475_good_list = [i475_good_list i];
            end
                

        end
    end

    
end

%%

i472_good_list = [];

for i = i472_expts
    if expt(i).TCs_extracted == 1 & ~isempty(expt(i).pass_run) & ~isempty(expt(i).dirtuning) 

        irun = 1;
        nrun = 1;
        run_str = ['runs']; 
        run_str = [run_str '-' expt(i).runs(irun,:)];

        LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        if strcmp(expt(i).folder,'lindsey')
            data_base = LG_base;
        elseif strcmp(expt(i).folder,'camaron')
            data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
        end

        mouse = expt(i).mouse;
        date = expt(i).date;
    
        load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
        if any(celleqel2mat_padded(adapt_input.tDoFeedbackMotion(21:end)))
            continue
        else
            ntrials = length(adapt_input.tDoFeedbackMotion);
            [s b] = selectCalc(adapt_input,21:ntrials);
            if s>=0.8 & abs(b)<0.1
                i472_good_list = [i472_good_list i];
            end
                

        end
    end

    
end

%%
cd('Z:\All_Staff\home\camaron\Analysis\2P')
save('good_expt_list.mat', 'i475_expts', 'i475_good_list', 'i472_expts', 'i472_good_list', 'i1402_expts', 'i1402_good_list', 'i1403_expts', 'i1403_good_list')

%%

function good_list = genExptListNoPassive(expts)
    dataset = 'oriAdapt_V1_cam';
    eval(dataset);
    good_list = [];
    % performance = [];
    
    for i = expts
        if expt(i).TCs_extracted == 1 & ~isempty(expt(i).dirtuning) 
    
            irun = 1;
            nrun = 1;
            run_str = ['runs']; 
            run_str = [run_str '-' expt(i).runs(irun,:)];
    
            LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
            if strcmp(expt(i).folder,'lindsey')
                data_base = LG_base;
            elseif strcmp(expt(i).folder,'camaron')
                data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
            end
    
            mouse = expt(i).mouse;
            date = expt(i).date;
        
            load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
            if any(celleqel2mat_padded(adapt_input.tDoFeedbackMotion(21:end)))
                continue
            else
                ntrials = length(adapt_input.tDoFeedbackMotion);
                [s b] = selectCalc(adapt_input,21:ntrials);
                if s>=0.8 & abs(b)<0.1
                    good_list = [good_list i];
    
                end
                    
    
            end
        end
    
        
    end
end
