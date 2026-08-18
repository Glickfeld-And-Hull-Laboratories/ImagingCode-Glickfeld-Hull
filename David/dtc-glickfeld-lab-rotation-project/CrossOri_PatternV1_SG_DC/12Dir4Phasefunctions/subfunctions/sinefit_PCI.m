function [b_hat, amp_hat, per_hat,pha_hat,sse,R_square] = sinefit(theta,data)
%% rearrange the input data file
theta=theta(:);
thetas = unique(theta);
ntheta = length(thetas);
if iscell(data)
    meandata = cellfun(@mean,data);
    
    data_new = cell2mat(data);
    theta_new = [];
    for i =1:numel(data)
    theta_new = [theta_new;repmat(theta(i,1),numel(data{i,1}),1)];
    end 
    [b,i]=sort(meandata,'descend');
    maxtheta = theta(i(1));
    theta=theta_new;
    data = data_new;
    
elseif length(theta) == ntheta
    meandata = data;
    [b,i]=sort(meandata,'descend');
    maxtheta = theta(i(1));
else
    meandata = zeros(1,ntheta);
    for itheta = 1:ntheta
        ind = find(theta == thetas(itheta));
        meandata(1,itheta) = mean(data(:,ind),2);
    end
    [b,i]=sort(meandata,'descend');
     maxtheta = thetas(i(1));    
end
    data=data(:);    




% make initial guesses, find the maxdirection tuning in the mean data to
% make u1_guess, incase similar too peaks, try another initial guess with a
% pi apart

amp_guess=(b(1)-b(end))/2;
b_guess=mean(b);
per_guess=2*pi;
ff= fft(meandata);
pha_guess= angle(ff(2));

if b_guess>0
    b_lb = b_guess.*.5;
    b_ub = b_guess.*2;
else
    b_ub = b_guess.*.5;
    b_lb = b_guess.*2;
end
amp_lb = 0;
amp_ub = amp_guess.*2;
per_lb = per_guess;
per_ub = per_guess;
pha_lb = -4.*pi;
pha_ub = 4.*pi;

% use fminsearch

options = optimset('MaxFunEvals',inf,'MaxIter',100000);
[out,fval,success] = fminsearchbnd(@sine_sse, [b_guess, amp_guess,per_guess,pha_guess],[b_lb,amp_lb, per_lb, pha_lb],[b_ub,amp_ub, per_ub, pha_ub], options);
if success == 1
    b_hat=out(1);
    amp_hat=out(2);
    per_hat=out(3);
    pha_hat=out(4);
    sse=fval;
    sse_tot = sum((data - mean(data)).^2);
    R_square = 1-(sse/sse_tot);
else
    b_hat=nan;
    amp_hat=nan;
    per_hat=nan;
    pha_hat=nan;
    sse=nan;
    sse_tot = nan;
    R_square = nan;
end


    % define the nested subfunction
    function miaosse = sine_sse(in)
        % pull out the slope and intercept
        b_tmp = in(1);
        amp_tmp = in(2);
        per_tmp = in(3);
        pha_tmp = in(4);
       
        y_fit = b_tmp+amp_tmp.*(sin(2*pi*theta./per_tmp + 2.*pi/pha_tmp));
        
        residuals = data - y_fit;
        miaosse = sum(residuals.^2);
    end
end


% y=vonmises(b,k1,k2,u1,theta);
% plot(theta,y)
% set(gca,'XTick',[0:pi/6:2*pi],'XTicklabel',directiondeg) 
