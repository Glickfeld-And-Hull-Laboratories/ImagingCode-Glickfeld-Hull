function order_of_days = generate_real_dates(years_array)
    months = ['01','02','03','04','05','06','07','08','09','10','11','12'];
    order_of_days = [];
    for y = years_array
        y = num2str(y);
        year = y(end-1:end);
        for m = 1:2:length(months)
            month = [months(m),months(m+1)];
            switch month
                case '01'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
                case '02'
                    if mod(str2num(year),4) == 0
                        days = [str2num([year,month,'01']) : str2num([year,month,'29'])];
                    else
                        days = [str2num([year,month,'01']) : str2num([year,month,'28'])];
                    end
                case '03'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
                case '04'; days = [str2num([year,month,'01']) : str2num([year,month,'30'])];
                case '05'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
                case '06'; days = [str2num([year,month,'01']) : str2num([year,month,'30'])];
                case '07'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
                case '08'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
                case '09'; days = [str2num([year,month,'01']) : str2num([year,month,'30'])];
                case '10'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
                case '11'; days = [str2num([year,month,'01']) : str2num([year,month,'30'])];
                case '12'; days = [str2num([year,month,'01']) : str2num([year,month,'31'])];
            end
            order_of_days = [order_of_days,days];
        end
    end
end