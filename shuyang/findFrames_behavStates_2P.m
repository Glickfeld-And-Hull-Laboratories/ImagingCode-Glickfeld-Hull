%find frames for running, stationary, back and forth, 
%used in behav_analysis_movingDots_2P

function[frames,frames_stay_cell, frames_bf_cell, frames_run_cell, frames_move_cell] = findFrames_behavStates_2P(speed)

frames = 1: length(speed);

%% generate a cell for stay parts
frames_stay_cell = {};
frames_stay = frames(speed == 0);
bound_stay = find(diff(frames_stay)~=1);
for j = 1:length(bound_stay)
    if j == 1
        frames_stay_cell{j} = frames(1:bound_stay(1));
    else
        frames_stay_cell = cat(2, frames_stay_cell, frames_stay(bound_stay(j-1)+1):frames_stay(bound_stay(j)));
    end
end
% delete the parts that are too short
transient_stay = 15;
for k = 1: size(frames_stay_cell,2)
    if length(frames_stay_cell{k}) <= transient_stay
        frames_stay_cell{k} = [];
    end
end
empties = find(cellfun(@isempty,frames_stay_cell));
frames_stay_cell(empties) = [];


%% find frames for move
frames_move = frames(speed ~= 0);
diff_frames_move = diff(frames_move);
gap = 6; % if the animal is still for less than 200ms during running, the running before and after the short still should still be counted as one part
bound_move = find(diff_frames_move>gap);
% put continuous frames together, generate a cell for all running parts
frames_move_cell ={};
for j = 1:length(bound_move)
    if j == 1
        frames_move_cell{j} = frames_move(1):frames_move(bound_move(j));
    else
        frames_move_cell{j} = frames_move(bound_move(j-1)+1):frames_move(bound_move(j));
    end
end

%delete the parts that are too short
framesMinDura = 6;
for k = 1: size(frames_move_cell,2)
    if length(frames_move_cell{k}) <= framesMinDura
        frames_move_cell{k} = [];
    end
end
empties = find(cellfun(@isempty,frames_move_cell));
frames_move_cell(empties) = [];

%% find back and forth parts (frames_bf_cell), and run parts
frames_bf_cell = {};
frames_run_cell = {};

for m = 1: size(frames_move_cell,2)
    bf = find(max(speed(frames_move_cell{m})) <= 31);
    if ~isempty(bf) &&  min(speed(frames_move_cell{m})) < 0
        frames_bf_cell = cat(2, frames_bf_cell, frames_move_cell{m});
    else
        frames_run_cell = cat(2, frames_run_cell, frames_move_cell{m});
    end
end

end





