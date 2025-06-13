function [s2d, p2d] = sin2d(x,ls,d,phi)
% function [s2d, p2d] = sin2d(x,ls,d,phi)
% x	[x y] values
% ls	wavelength
% d	phase % radians
% phi	ori  in degrees

% Returns the phase of the response (p2d) and the sin of the response (s2d)

k = [cos(pi*phi/180) sin(pi*phi/180)];
v = x;
kp = [cos((phi+90)*pi/180) sin((phi+90)*pi/180)];
%size(k)
par = v*k';
%size(par)
perp = v*kp';
%size(perp)
p2d = (2*pi*par/ls) + d;
s2d =  (cos((2*pi*par/ls) + d));
