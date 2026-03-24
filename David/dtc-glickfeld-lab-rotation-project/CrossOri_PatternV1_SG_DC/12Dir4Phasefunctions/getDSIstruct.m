% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)

function [DSIstruct] = getDSIstruct(avg_resp_dir)
    
    nCells  = size(avg_resp_dir,1);
    nDir    = size(avg_resp_dir,2);
    for iCell = 1:nCells
        [max_val max_ind]   = max(avg_resp_dir(iCell,:,1,1,1));
        null_ind            = max_ind+(nDir./2);
        null_ind(find(null_ind>nDir)) = null_ind(find(null_ind>nDir))-nDir;
        min_val     = avg_resp_dir(iCell,null_ind,1,1,1);
        if min_val < 0; min_val = 0; end
        DSI(iCell)          = (max_val-min_val)./(max_val+min_val);
        DSI_maxInd(iCell)   = max_ind; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Global DSI & Global OSI
    angs = 0:(360/nDir):(360 - 360/nDir);
    
    for j = 1:nCells
        amps = avg_resp_dir(j,:,1,1,1); %our responses   
        
        g_dsi(j) = sqrt( sum(sin(1*angs*pi/180).*amps).^2 + sum(cos(1*angs*pi/180).*amps).^2)/sum(amps);
        g_osi(j) = sqrt(sum(sin(2*angs*pi/180).*amps).^2 + sum(cos(2*angs*pi/180).*amps).^2)/sum(amps);
        
        xm_dsi = (sum(amps.*cos(deg2rad(1*angs)))/sum(amps)); %mean of the response, x
        ym_dsi = (sum(amps.*sin(deg2rad(1*angs)))/sum(amps)); %mean of the response, y
        xm_osi = (sum(amps.*cos(deg2rad(2*angs)))/sum(amps));
        ym_osi = (sum(amps.*sin(deg2rad(2*angs)))/sum(amps));
        
        ang_dsi(j) = rad2deg(atan2(ym_dsi, xm_dsi)); %preferred direction, in degrees
        ang_osi(j) = (rad2deg(atan2(ym_osi, xm_osi)));
        ang_osi(j) = ang_osi(j)/2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    DSIstruct.DSI           = DSI;  % DSI value, ranging 0 to 1
    DSIstruct.DS_ind        = find(DSI>0.5);    % index of cells meeting criteria DSI>0.5
    DSIstruct.prefDir       = DSI_maxInd;   % index of preferred direction, ranging 1 through nDir
    DSIstruct.gDSI          = g_dsi;    %global DSI value, from 0 to 1
    DSIstruct.gDSI_prefDir  = ang_dsi;  %preferred direction, calculated from gDSI
    DSIstruct.gOSI          = g_osi;    %global DSI value, from 0 to 1
    DSIstruct.gOSI_prefDir  = ang_osi;  %preferred direction, calculated from gDSI

end