
function [osi,ang] = compdsi(angs,amps)
% function [osi,ang] = composi(angs,amps)
% function dsi = compdsi(angs,amps)
% angs is in degrees, amps in whatever units

osi = sqrt( sum(sin(angs*pi/180).*amps).^2 + sum(cos(angs*pi/180).*amps).^2)/sum(amps);
%osi = 1- abs(sum(amps.*exp(2*i*(angs)))/(sum(amps)));
%num = abs(sum(amps.*exp(2*i*(angs))))
%denom = sum(amps)

xm = (sum(amps.*cos(deg2rad(angs)))/sum(amps));
ym = (sum(amps.*sin(deg2rad(angs)))/sum(amps));

r = sqrt(xm^2 + ym^2);

ang = (rad2deg(atan(ym/xm)));

if xm<0
  ang = 180+ang;
elseif ym<0
  ang = 360+ang;
end

ang;



