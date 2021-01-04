clear;
folder = 'Y:\home\nathan\DATA\S-probe\Awake\behavior\';
sessionID =  '1045-191210';
filename = dir([folder 'data-i' '*' sessionID  '*' ]);
file = [folder, filename.name];
load(file);
countTime = cell2mat(input.counterTimesUs);
countValue = cell2mat(input.counterValues); % this variable tells you which pulse it is from master8
quadTime = cell2mat(input.quadratureTimesUs); % these 2 variables are for speed
quadValue = cell2mat(input.quadratureValues);
% 1. calculate speed:
speed = [];
nStart = find(countValue == 1, 1, 'last'); % find the real start of the experiment (get rid of the press button twice shit)
for n = nStart:length(countTime)            % you can do the same thing as this for loop when counting spikes
    if countValue(n) == 1                    
        this_frame_quad_inx = find(quadTime > countTime (n) - 100*1000 & quadTime <= countTime(n));
    else
        this_frame_quad_inx = find(quadTime > countTime(n-1) & quadTime <= countTime(n));
    end
    this_frame_loc = quadValue(this_frame_quad_inx);
    this_frame_distance = this_frame_loc(end)-this_frame_loc(1);
    this_frame_time = quadTime(this_frame_quad_inx);
    this_frame_timelength = this_frame_time(end)-this_frame_time(1);
    this_frame_speed = (this_frame_distance)*1000000/(this_frame_timelength);
    % 1000000: us to s
    speed = [speed this_frame_speed]; % using cat is not the most efficient way, it will take a few minutes to run but not too bad
end

airpuffon = input.cTactileStimTurnedOn;
airpuffon1 = airpuffon(~cellfun('isempty',airpuffon)); %it seems the airpuff sometimes is skipped, and the cell is empty, through away these if any
airpuffon1 = double(cell2mat(airpuffon1)); % this tells you which pulse the airpuff starts.
lenairpuff = double(input.tactileStimulusDurationUs);
lenairpuff = lenairpuff/1000000; % length of airpuff in seconds

% get behavioral states
frm_maxGap = 3; %300ms in between
[frames,frames_stay_cell, frames_bf_cell, frames_run_cell, frames_move_cell] = findFrames_behavStates(speed,frm_maxGap);
frames_stay = cell2mat(frames_stay_cell);
frames_run = cell2mat(frames_run_cell);
behav_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\behavior';
save([behav_dest '\' sessionID '_behavAnalysis.mat' ],...
    'speed','frames','frm_maxGap','frames_stay_cell','frames_bf_cell',...
    'frames_run_cell','frames_move_cell','frames_stay','frames_run');

%%
sessionID = '1046-191216';
behav_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\behavior';
behav_struct = load([behav_dest '\' sessionID '_behavAnalysis.mat' ]);
frames_run_cell = behav_struct.frames_run_cell;
speed = behav_struct.speed;
[frames_befo_run_cell,frames_aft_run_cell,frames_runTrigger_mat,...
    frames_runoff_include,frames_runoff_mat,frames_run_mat] = findFrames_runWindows (speed,frames_run_cell);
save([behav_dest '\' sessionID '_behavAnalysis.mat' ],...
    'frames_befo_run_cell','frames_aft_run_cell','frames_runTrigger_mat',...
    'frames_runoff_mat','frames_run_mat', 'frames_runoff_include','-append');



