function y2 = interp_frame(tc_tot, end2)
% tc is ncell by nFrame
y2= [];
for id = 1: size(tc_tot,2)  % for each session
    tc = tc_tot{id};
    if isnan(tc)
        continue
    else
        end1 = size(tc,2); %nFrames
        y2_temp = zeros(size(tc,1),end2); %ncells by frames
        for i = 1:size(tc,1) %nCells
            
            y2_temp(i,:) = interp1(1:1:end1, tc(i,:), 1:(end1 - 1)/(end2-1):end1);
        end
        y2 = [y2; y2_temp]; %ncells by frames
    end
end
end