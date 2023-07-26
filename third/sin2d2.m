function [s2d, p2d] = sin2d2(xs,ls,d,phi)
% function [s2d, p2d] = sin2d(xs,ls,d,phi)
% xs	x and y values  % assume square
% ls	wavelength
% d	phase (radians)
% phi	ori  (degrees)

k = [cos(pi*phi/180) sin(pi*phi/180)];
kp = [cos((phi+90)*pi/180) sin((phi+90)*pi/180)];
xvec = xs;
yvec = xs;

p2d  = zeros(length(xvec),length(xvec));
s2d  = zeros(length(xvec),length(xvec));

for j=1:length(xvec)
  for l=1:length(yvec)
    	v = [xvec(j) yvec(l)];
	par = v*k';
	perp = v*kp';
	p2d(j,l) = (2*pi*par/ls) + d;
	s2d(j,l) =  (cos((2*pi*par/ls) + d));
  end
end

