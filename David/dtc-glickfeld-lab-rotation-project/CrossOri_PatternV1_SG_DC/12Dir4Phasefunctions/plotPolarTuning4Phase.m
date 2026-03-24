% Input arguments from function *getAlignedGratPlaidTuning.m*
% where both grating and plaid tuning curves are aligned to "preferred"
% grating

function plotPolarTuning4Phase(gratAligned, plaidAligned, iCell, phases, plot_type)
    
    if nargin < 4   % If only the first three inputs, assume all phases are plotted
        plot_type = 'default';
    end

    x       = -150:30:180;
    x_rad   = deg2rad(x);

    colors = getColors;

    if strcmp(plot_type, 'separate') %this just allows us to either plot the responses separately or together on the same figure
        figure(555+iCell)
        plaid_max = max(plaidAligned(iCell,:,:), [], 'all');
        grat_max = max(gratAligned(iCell,:), [], 'all');
        max_r = max(plaid_max, grat_max);
        for im = 1:4
            subplot(1,4, im)
            polarplot([x_rad x_rad(1)], [plaidAligned(iCell,:,im) plaidAligned(iCell,1,im)],'Color',colors(im,:), 'LineWidth',3.5)
            hold on
            polarplot([x_rad x_rad(1)], [gratAligned(iCell,:) gratAligned(iCell,1)],'k', 'LineWidth',3.5)
            ax = gca;
            ax.RLim = [0 max_r];
        end
    elseif strcmp(plot_type,'together')
        figure(555+iCell)
        for phase = phases
            polarplot([x_rad x_rad(1)], [plaidAligned(iCell,:,phase) plaidAligned(iCell,1,phase)],'Color',colors(phase,:), 'LineWidth',4.5)
            hold on
        end
        polarplot([x_rad x_rad(1)], [gratAligned(iCell,:) gratAligned(iCell,1)],'k', 'LineWidth',4.5)
    elseif strcmp(plot_type,'default')
        for im = 1:4
        polarplot([x_rad x_rad(1)], [plaidAligned(iCell,:,im) plaidAligned(iCell,1,im)],'Color',colors(im,:))
        hold on
        end
        polarplot([x_rad x_rad(1)], [gratAligned(iCell,:) gratAligned(iCell,1)],'k', 'LineWidth',2) 
    end
end