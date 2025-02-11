function [reporter, reporter1] = xcorrBusiness(unit1, unit2, bin, xmin, xmax)

numbins = (xmax - xmin)/bin;
counts = zeros(numbins,1);
unit2 = unit2*1000; % convert to ms
unit2 = floor(unit2);
reporter1 = unit2;
unit2train = zeros(unit2(end),1);
for i = 1:length(unit2)
    unit2train(unit2(i)) = 1;
    reporter = unit2train;
end
k = 0;
for i = 21:(length(unit1)-21)
  TS = floor(unit1(i)*1000);
  reporter1(i) = TS;
  for j = 1:(numbins-1)
      
      if unit2train(TS(i),1) %((TS + xmin)+(i-1)) == 1
          counts(j,1) = counts(i)+1;
          k = k+1;
      
      end
  end
  reporter = k;
end
bar(counts)
          
      
         
    
