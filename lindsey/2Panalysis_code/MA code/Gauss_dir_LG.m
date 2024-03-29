%VB Gabor fitting, modified to fit sum of gaussians
function h = Gauss_dir_LG(pars,angle)

 r0 = pars(1); %offset amp
 rp = pars(2); %pref amp
 rp2 = pars(3); %null amp
 theta = angle; % rad
 thetap = pars(4); % rad
 sigma = pars(5); %rad
 if thetap<pi
    h = r0 + rp.*exp(-1*((theta-thetap).^2)./(2.*(sigma^2)))+ rp2.*exp(-1*((theta-thetap-pi).^2)./(2.*(sigma^2)));
 else
    h = r0 + rp.*exp(-1*((theta-thetap).^2)./(2.*(sigma^2)))+ rp2.*exp(-1*((theta-thetap+pi).^2)./(2.*(sigma^2)));
 end

% from VonMises:
% if thetap<pi
%     h = r0 + rp * exp( 2* (cos(theta-thetap) -1 )) + rp2 * exp(2 * (cos(theta-thetap-pi) -1 )) ;
% else
%     h = r0 + rp * exp( 2* (cos(theta-thetap) -1 )) + rp2 * exp(2 * (cos(theta-thetap+pi) -1 )) ;
% end

