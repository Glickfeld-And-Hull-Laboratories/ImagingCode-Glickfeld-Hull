function x = greens(num,scale);
% create a color map of blue: pale to dark
x = zeros(num,3);
x(:,2) = 1;
if nargin<2
    scale = 'linear';
end
if strcmp(scale,'log')
    x(:,1) = 1-logspace(log(1/num),log(1-(1/num)),num);
    x(:,3) = 1-logspace(log(1/num),log(1-(1/num)),num);
elseif strcmp(scale,'linear')
    x(:,1) = 1-(1/num):1/num:0;
    x(:,3) = 1-(1/num):1/num:0;
end