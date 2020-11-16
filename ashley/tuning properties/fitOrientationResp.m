function [fit,R2] = fitOrientationResp(avgResp,orientations)
    
    ori_rad = deg2rad(orientations);
    theta_smooth = 0:1:179;

    [b, k, R,u,~,R2] = ...
            miaovonmisesfit_ori(ori_rad,avgResp);
    fit = b+R.*exp(k.*(cos(2.*(deg2rad(theta_smooth)-u))-1));
    
end