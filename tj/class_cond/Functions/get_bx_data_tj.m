function b_data = get_bx_data_tj(bdata_source, session);
%if statement to distinguish between 1000s, 900s, and 00s and find bx file

%get animal number and session date
img_loc = strfind(session, 'img');
animal_num = session(img_loc+3:img_loc+6);
session_date = session(1:6);
session_time = session(img_loc+8:img_loc+11)

%identify behavior file

bfile = fullfile(bdata_source, ['data-i' animal_num  '-' session_date '-' session_time '.mat' ]);

%load behavior file
assert(length(bfile)) = 1;
b_data = load(bfile);
b_data = b_data.input;
end

