% function designed to determine that state of the pockel cell by using the changes in whole field fluorescence
function pockel_tc = find_pockel_tc(data, tc_switch)
%input - data - should be the full 3D matrix of the 2P imaging data
%tc_switch - should =0 for full data movie and 1 for just the TCs  size = frames,neurons
%output - pockel_tc - a boolean timecourse where frames with laser power =1    and frames where pockel cell 
% reduces laser power = 0

%calculate grand avg
if tc_switch == 1
    data_grand_avg = mean(data,2)';
else
    if numel(size(data)) > 3
        data=squeeze(data);
    end
    data_grand_avg = squeeze(mean(mean(data,2),1))'; % average across vertical and horizontal pixels, each frame is a number
end

%calculate some useful variables
data_round = round(data_grand_avg);
data_mode = mode(data_round(1:round(length(data_grand_avg)/3))); % most frequent values in data_avg
data_deriv = diff(data_grand_avg);
data_deriv_std = std(data_deriv);

%use the derivative to find transitions between power up and power down
[~, laser_on] = findpeaks(data_deriv, 'MinPeakHeight', data_deriv_std*5);
[~, laser_off] = findpeaks(data_deriv*-1, 'MinPeakHeight', data_deriv_std*5);

%if laser starts powered up then make all frames until the end of the first trial = 1
pockel_tc = zeros(1,length(data_grand_avg));
if laser_on(1) > laser_off(1) %imaging starts with laser power down
    pockel_tc(1,1:(laser_off(1)-2)) = 1;
    %laser_off = laser_off(2:end);
else
    %verify that the laser starts in the off condition
    %assert(mean(data_grand_avg(3:10)) < data_mode*1.25);
end

%use transition frame indeces to calculate laser power (pockel) timecourse

for laser_trans_ind = laser_on
    if laser_off(end) > laser_trans_ind
        pockel_tc(1,laser_trans_ind+2:laser_off(find(laser_off>laser_trans_ind, 1, 'first'))-2) = 1; %remove in two frames of buffer on either side
    else
        pockel_tc(1,laser_trans_ind+2:end) = 1;
    end
end
%pockel_tc(1) = 0% added this because first frame is a transition, satisfies next assertion
%verify that all laser_on frames are indeed powered
%assert( min(data_grand_avg(find(pockel_tc))) > data_mode );

return


