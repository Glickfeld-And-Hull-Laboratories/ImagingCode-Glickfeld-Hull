function [this_mouse, this_date] = select_mouse_info(sess_info);

%get animal number 
img_loc = strfind(sess_info, 'img');
this_mouse = sess_info(img_loc+3:end);
if length(this_mouse) <3
    this_mouse = strcat('9', this_mouse);  %for the behavior analysis and some others the 900s omit the first digit in the animal number so img955 would be img55
end

%get session date
this_date = sess_info(1:6);

return