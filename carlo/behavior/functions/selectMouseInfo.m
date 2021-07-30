function [thisMouse, thisDate] = selectMouseInfo(sess_info, dateIdx) %(info)

%get animal number 
img_loc = strfind(sess_info, 'img');
thisMouse = sess_info(img_loc+3:end);
if length(thisMouse) <3
    thisMouse = strcat('9', thisMouse);  %for the behavior analysis and some others the 900s omit the first digit in the animal number so img955 would be img55
end

%get session date

thisDate = sess_info(1:dateIdx);


% this_mouse = num2str(input.subjectNum);
% this_date = input.saveTime;

return