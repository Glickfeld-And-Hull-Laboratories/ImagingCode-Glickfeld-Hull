k = 1;
cofire = 0;
for i = 1:500
timeInterest1(i) = unit497(i) - .005;
timeInterest2(i) = unit497(i) + .005;
cotimes = unit473((unit502 < timeInterest2(i)) & (unit502 > timeInterest1(i)));
for j = 1:length(cotimes)
    cofire(k) = cotimes(j);
    k = k + 1;
end

end