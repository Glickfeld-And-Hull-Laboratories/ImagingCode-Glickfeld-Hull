plane_number = 1; %1 => 1st plane (or volume) imaged that day

PWD = 'G:\users\lindsey\analysisLG\active mice\AC41\110728\';
PWD_EEG = 'G:\users\lindsey\dataLG\2010\110210\LG25\';

file_EEG_mat = strvcat(...
    'eegacq_mouse_training2011-2-9_16-11-40_rec_110209_CFJ_run1.mat', ...
    'eegacq_mouse_training2011-2-9_16-49-43_rec_110209_CFJ_run2.mat', ...
    'eegacq_mouse_training2011-2-9_17-19-20_rec_110209_CFJ_run3.mat', ...
    'eegacq_mouse_training2011-2-9_17-55-15_rec_110209_CFJ_run4.mat', ...
    'eegacq_mouse_training2011-2-9_18-24-37_rec_110209_CFJ_run5.mat');

file_mat = strvcat(...
    '110728_AC41_run1.tif', ...
    '110728_AC41_run2.tif', ...
    '110728_AC41_run3.tif');

%file_mat = ['100502_u12_run3023__ch2_fix_reg.tif'];
file_ch2 = 'AVG_110210_LG25_run1.tif';
file_mask = 'AVG_110210_LG25_run1.tif';

DIR_rand0 = strvcat('G:\users\lindsey\dataLG\vis_stim_protocols\110728\AC41_run2\TF50_SF50_ori50_Con50_Nfrpstim240_Ret1_GabSz20_28Jul2011\');
                
Nprotocols = size(DIR_rand0,1);

randvec = [1 1; ...
           1 1]; %directory #, randnum (within that directory)
           

RandStim_info = [];

%%%Repeat entry of these parameters for each stimulus protocol presented (i.e. Nprotocols): 
count_prot = 1;
RandStim_info.prot(count_prot).userun0 = [1];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [repmat([2*ones(1,3) 8*ones(1,3)],1, 8) 2 8 0];
RandStim_info.prot(count_prot).SFvec = [repmat([.06 .12 .24], 1, 16) 0 0 0];
RandStim_info.prot(count_prot).pos = [ones(1,51)];
RandStim_info.prot(count_prot).StimXpos = [10*ones(1,50) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-10*ones(1,50) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 50; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [0*ones(1,6) 45*ones(1,6) 90*ones(1,6) 135*ones(1,6) 180*ones(1,6) 225*ones(1,6) 270*ones(1,6) 315*ones(1,6) 0 0 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,50) 0];
RandStim_info.prot(count_prot).Fs1 = 4; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 4; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 19;
RandStim_info.prot(count_prot).Non1 = 19;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);


%%% in case of second protocol on the same day, fill in the following: 
if Nprotocols > 1
    count_prot = 2;
RandStim_info.prot(count_prot).userun0 = [5];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [2 8 2*ones(1,8) 8*ones(1,8) 0];
RandStim_info.prot(count_prot).SFvec = [0.02 0.02 0.32*ones(1,16) 0];
RandStim_info.prot(count_prot).pos = [ones(1,19)];
RandStim_info.prot(count_prot).StimXpos = [-5*ones(1,18) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-15*ones(1,18) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 =19; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [0 0 0 45 90 135 180 225 270 315 0 45 90 135 180 225 270 315 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,18) 0];
RandStim_info.prot(count_prot).Fs1 = 4; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 4; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 18;
RandStim_info.prot(count_prot).Non1 = 18;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end

if Nprotocols > 2
    count_prot = 3;
RandStim_info.prot(count_prot).userun0 = [6];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [8*ones(1,5) 0];
RandStim_info.prot(count_prot).SFvec = [0.02 0.04 0.08 0.16 0.32 0];
RandStim_info.prot(count_prot).pos = [ones(1,6)];
RandStim_info.prot(count_prot).StimXpos = [-10*ones(1,5) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [0*ones(1,5) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 6; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [0*ones(1,5) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,5) 0];
RandStim_info.prot(count_prot).Fs1 = 4; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 4; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 18;
RandStim_info.prot(count_prot).Non1 = 18;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end
%%%%% next, enter all the info above so it's easy to access later if needed
RandStim_info.general.PWD = PWD;
RandStim_info.general.PWD_EEG = PWD_EEG;
RandStim_info.general.file_EEG_mat = file_EEG_mat;
RandStim_info.general.file_mat = file_mat;
RandStim_info.general.file_ch2 = file_ch2;
RandStim_info.general.file_mask = file_mask;
RandStim_info.general.DIR_rand0 = DIR_rand0;
RandStim_info.general.Nprotocols = size(DIR_rand0,1);
RandStim_info.general.randvec = randvec;
RandStim_info.general.plane_number = plane_number; 
save([PWD,'RandStim_info_plane',num2str(plane_number),'.mat'],'RandStim_info');

