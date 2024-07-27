

% for establishing the following info from recordings during experiments
    % probe/headstage combination
    % genotype
    % experiment list
    % spike sorted file name
    
%%
clear params
save_path = 'Z:\home\Naomi\Electrophysiology\Kilosort_Analysis\'
probeLog_path = 'Z:\home\Naomi\Electrophysiology\New Analysis\probeLogNB.xlsx';

expdate=input('Date?')
params.info.recDate = expdate;
params.info.probeType = input('Probe Type? (A32=3, NP=4)');
% types =
    % (1) H4/optrode
    % (2) A32/optrode
    % (3) A32
    % (4) Neuropixel
params.info.headstage = input('Headstage? (BR=2, NP=4)');
% types = 
    % (1) old blackrock (pre 12/2021)
    % (2) new blackrock same side omnetics  
    % (3) new blackrock opp side omnetics (uncommon)
    % (4) NP Headstage
display('All are characters')
params.info.Genotype = input('Genotype? ');
params.info.ID_Num = input('ID? ')
params.info.PDay= input('Age? ')
params.info.ExpSize= input('Field? ')
params.info.sex= input('Sex? ')
params.info.cog= input("Awake or Anest?")
params.info.expt_list = input('Variables? (Use Brackets)')

%%
% function wrapExptParams(params,save_path,probeLog_path)
probeList = {'H4-optrode';'A32-optrode';'A32';'Neuropixel';};
hsList = {'old Blackrock';'new Blackrock (same side)';'old Blackrock (opposite side)';'NP Headstage'};

params.info.probeType = probeList{params.info.probeType};
params.info.headstage = hsList{params.info.headstage};
%%
filenum=input('Which file for spiketimes? (Use hangers)')
params.info.spikeFiles = (['datafile',filenum]);
%%
save(fullfile(save_path,params.info.recDate,[params.info.recDate,'_',params.info.spikeFiles,'_exptInfo.mat']))

disp(['wrapped experiment info for ',params.info.recDate, ' :)']);
% end
