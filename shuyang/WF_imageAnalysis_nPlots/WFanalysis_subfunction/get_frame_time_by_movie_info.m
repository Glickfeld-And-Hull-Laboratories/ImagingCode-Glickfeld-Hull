function [ frame_times ] = get_frame_time_by_movie_info( info )
% get the frame timing information from movie information
% This function has tons of assumptions and should be used with extreme
% caution

INFO_ID = 51123; % this is probabaly camera specific!!!
frame_times = zeros(length(info),1);
for i=1:length(info)
    for j=1:length(info(i).UnknownTags)  % look for unknown tag
        if(info(i).UnknownTags(j).ID ==INFO_ID)
            %extract the timing information
            str = info(i).UnknownTags(j).Value;
            pre_str = '"ElapsedTime-ms":';
            inx = strfind(str, pre_str);
            % find next comma
            reg_float = '[-+]?(\d*[.])?\d+';
            reg_exp = [pre_str reg_float];
            full_time_str = regexp(str,reg_exp ,'match');
            time_str = regexp(full_time_str{1},reg_float ,'match');
            frame_times(i) = str2num(time_str{1});
            break;
        end
        assert(j < length(info(i).UnknownTags), 'did not find Value');
    end
end
end

