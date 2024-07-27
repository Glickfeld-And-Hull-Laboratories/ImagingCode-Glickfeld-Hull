function [ys] = smooth(y,sd)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if size(y,2)==1
    y=y';
end
[Nr,Nc]=size(y);
%Create Gaussian
sd2=sd*sd;
Ng=4*sd;
if Ng<1
    Ng=1;
end
xg=[-Ng:Ng];
g=exp(-0.5*(xg.^2)/sd2);
g=g/sum(g);
ys = zeros(size(y));
for j=1:Nr
    z=convn(y(j,:),g,'same');
    ys(j,:)=z;
end

end

