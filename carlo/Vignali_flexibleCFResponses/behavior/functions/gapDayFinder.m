function gapInx = gapDayFinder(days);
%find the gaps where training days were not consecutive. For a list of days this function convets them to numbers and returns the values between the
%indeces of non-consecutive days. e.g. 2.5 if days 2 and 3 were not back to back. 

%establish order of days in a year. 
%{
order_of_days = cat(2, [170101:170131],[170201:170228],[170301:170331], [170401:170430], [170501:170531], [170601:170630], ...
    [170701:170731], [170801:170831], [170901:170930], [171001,171031], [171101:171130], [171201:171231],      [180101:180131],[180201:180228],[180301:180331], [180401:180430], [180501:180531], [180601:180630], ...
    [180701:180731], [180801:180831], [180901:180930], [181001,181031], [181101:181130], [181201:181231], ...
    [190101:190131], [190201:190228], [190301:190331], [190401:190430], [190501:190531], [190601:190630], [190701:190731], [190801:190831], [190901:190930], [191001:191031], [191101:191130], [191201:191231]);
%}
order_of_days = generate_real_dates([2019,2020,2021,2022,2023]);

gapInx = [];

for ii = 1:length(days)-1
    %convert days info to numbers. Get index of analyzed day
    this_day = str2double(days{ii}(1:6));
    next_day = str2double(days{ii+1}(1:6));
    this_day_inx = find(order_of_days==this_day);
    
    %ensure that the dates listed are included in the list of days
    assert(ismember(this_day, order_of_days));
    assert(ismember(next_day, order_of_days));
    
    %determine if current day and next day are consecutive
    if next_day == order_of_days(this_day_inx+1)
        continue
    else
        gapInx = [gapInx, ii+0.5];
    end
end