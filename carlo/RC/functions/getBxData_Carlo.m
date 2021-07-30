function b_data = get_bx_data_sj(bdata_source, session);
%if statement to distinguish between 1000s, 900s, and 00s and find bx file

%get animal number and session date
img_loc = strfind(session, 'img');
animal_num = session(img_loc+3:img_loc+6);
session_date = session(1:6);

%identify behavior file
if strlength(session) < 15
    bfile = dir([bdata_source 'data-i' animal_num '-' session_date '*' ]);
else
    bfile = dir([bdata_source 'data-i' animal_num '-' session_date session(end-4:end) '.mat']);
end

%load behavior file
behave_dest = [bdata_source bfile.name];
assert(length(bfile)) = 1;
b_data = load(behave_dest);
b_data = b_data.input;
end