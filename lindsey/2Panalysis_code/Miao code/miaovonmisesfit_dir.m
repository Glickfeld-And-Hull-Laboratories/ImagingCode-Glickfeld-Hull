function [b_hat,k1_hat,R1_hat,R2_hat,u1_hat,u2_hat,sse,R_square] = miaovonmisesfit_dir(theta,data)
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
    [b,i]=sort(meandata,'descend');
    maxtheta = thetas(i(1));
end
    data=data(:);    




% make initial guesses, find the maxdirection tuning in the mean data to
% make u1_guess, incase similar too peaks, try another initial guess with a
% pi apart


b_guess=b(end);
k1_guess=0.5;
R1_guess=b(1).*1.2;
u1_guess=maxtheta;% try the first peak 
orthogtheta = maxtheta+pi;
if orthogtheta>=(2.*pi)
    orthogtheta = orthogtheta- (2.*pi);
end
R2_guess = data(find(chop(thetas,4) == chop(orthogtheta,4)));

b_lb = 0;
b_ub = b(1);
k1_lb = .3466;
k1_ub = 30;
R1_lb = 0;
R1_ub = b(1);
R2_lb = 0;
R2_ub = R2_guess;
u1_lb = maxtheta-(pi/ntheta);
u1_ub = maxtheta+(pi/ntheta);

options = optimset('MaxFunEvals',inf,'MaxIter',100000);
[out,fval,success] = fminsearchbnd(@vonmises_sse, [b_guess, k1_guess,R1_guess,R2_guess,u1_guess],[b_lb,k1_lb, R1_lb, R2_lb, u1_lb],[b_ub,k1_ub, R1_ub, R2_ub, u1_ub],options);
if success == 1
    b_hat=out(1);
    k1_hat=out(2);
    R1_hat=out(3);
    R2_hat=out(4);
    u1_hat=out(5);
    u2_hat=mod(u1_hat+pi,2*pi);
    sse=fval;
    sse_tot = sum((data - mean(data)).^2);
    R_square = 1-(sse/sse_tot);
else
    b_hat=nan;
    k1_hat=nan;
    R1_hat=nan;
    R2_hat=nan;
    u1_hat=nan;
    u2_hat=nan;
    sse=nan;
    sse_tot = nan;
    R_square = nan;    
end


    % define the nested subfunction
    function miaosse = vonmises_sse(in)
        b_tmp = in(1);
        k_tmp = in(2);
        R1_tmp = in(3);
        R2_tmp = in(4);
        u1_tmp = in(5);

        y_fit = b_tmp+R1_tmp.*exp(k_tmp.*(cos(theta-u1_tmp)-1))+R2_tmp.*exp(k_tmp.*(cos(theta-u1_tmp-pi)-1));
        
        residuals = data - y_fit;
        miaosse = sum(residuals.^2);
    end
end


% y=vonmises(b,k1,k2,u1,theta);
% plot(theta,y)
% set(gca,'XTick',[0:pi/6:2*pi],'XTicklabel',directiondeg) 
