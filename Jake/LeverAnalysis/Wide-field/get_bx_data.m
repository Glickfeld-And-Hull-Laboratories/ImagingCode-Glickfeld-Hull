function b_data = get_bx_data(bx_source, day);
%if statement to distinguish between 1000s, 900s, and 00s and find bx file

%get animal number and session date
img_loc = strfind(day, 'img');
animal_num = day(img_loc+3:end);
if length(animal_num) <3
    animal_num = strcat('9', animal_num)  %for the behavior analysis and some others the 900s omit the first digit in the animal number so img955 would be img55
end
session_date = day(1:6);
img_loc
%identify behavior file
bfile = dir([bx_source 'data-i' '*' animal_num '-' session_date '*' ]);

%load behavior file
behave_dest = [bx_source bfile.name];
assert(length(bfile)) = 1;
b_data = load(behave_dest);
b_data = b_data.input;
end