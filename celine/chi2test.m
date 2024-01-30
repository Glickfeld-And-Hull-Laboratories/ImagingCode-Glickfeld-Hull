%?cs(data)
%Chi-Square Test =
%                       490
%Critical Value =
%          552.074742731856

% n = 5000
% k = 500
% n/k = 10

%Since 490 < 552, We accept the null hypothesis H_0 and conclude
%that we don't have enough evidence at level alpha = 0.05 to say that
%the observation are not uniform.

function chi2test(data, num_int)
	[m,n] = size(data);
	
	E = m/num_int; %E=10
	interval=zeros(num_int, 1);
	for i=1:m
		k = data(i,1);
		interval(fi(num_int,k),1) = interval(fi(num_int,k),1) +1;
	end
	x = 0;
	for j= 1:num_int
		x=x+ (interval(j,1)-E)^2/E;
	end
	x
	chi2inv(0.95,num_int-1)
	
function b = fi(i,n)
for j= 0:i
	if n > (1/i)*(i-j)
		b = i-j+1;
		break
	end
end


