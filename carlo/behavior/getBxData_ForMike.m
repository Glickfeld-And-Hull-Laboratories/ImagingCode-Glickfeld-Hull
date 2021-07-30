function bxData = getBxData_ForMike(bxSourceBase, day);
%if statement to distinguish between 1000s, 900s, and 00s and find bx file

%get animal number and session date
img_loc = strfind(day, 'img');
animal_num = day(img_loc+3:img_loc+6);
session_date = day(1:6);

%identify behavior file
if strlength(day) < 15
    bfile = dir([bxSourceBase 'data-i' animal_num '-' session_date '*' ]); 
else
    bfile = dir([bxSourceBase 'data-i' animal_num '-' session_date day(end-4:end) '.mat']);
end

%load behavior file
behave_dest = [bxSourceBase bfile.name]; %find path of bx data
assert(length(bfile)) = 1; %if more than one file is loaded, an error message is displayed
bxData = load(behave_dest); %satisfying the above condition; the dataset is loaded
bxData = bxData.input; %isolates input branch of struct (only branch, so idk why this is necessary - maybe to separate data during imagine)
end