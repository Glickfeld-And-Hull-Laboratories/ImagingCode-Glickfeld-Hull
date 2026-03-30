function [] = plot_phaseFit(ZpZcStruct, phaseModStruct, cell_indices)

phase_range = 0:359;
data_phase = [0 90 180 270];
PCI = phaseModStruct.PCI; %zpzc for every cell

colors = [
    212, 26, 75;  % pink        
    0, 174, 239;   % blue 
    237, 177, 32  % yellow
    126, 47, 142;  %  purple
    ] / 255;  % Normalize to [0,1]

if numel(cell_indices) > 1
    num_cells = numel(cell_indices);
    ncolumns = ceil(sqrt(num_cells));
    nrows = ceil(num_cells / ncolumns);
    
    figure
    for idx = 1:num_cells
        celln = cell_indices(idx);
        
        zpzc = PCI(:, celln);
        sine_fit = phaseModStruct.yfit(celln,:);
        
        subplot(nrows, ncolumns, idx)
        hold on
        for k = 1:4
            plot(data_phase(k), zpzc(k,:),'o','color', colors(k,:), 'MarkerFaceColor', colors(k,:));
        end
        plot(phase_range, squeeze(sine_fit), 'k-', 'LineWidth', 1.5)
        title(sprintf('Cell %d', celln));
        ylim([-7 7]);
        xlabel('phase')
        ylabel('Zp - Zc')

    end
else
    figure
    hold on
    for celln = cell_indices
        zpzc = PCI(:, celln);
        sine_fit = phaseModStruct.yfit(celln,:);
        
        for k = 1:4
            plot(data_phase(k), zpzc(k,:),'o','color', colors(k,:), 'MarkerFaceColor', colors(k,:));
        end
        plot(phase_range, sine_fit, 'k-', 'LineWidth', 1.5)
        title(sprintf('Cell %d', celln));
        ylim([-7 7]);
        xlabel('phase')
        ylabel('Zp - Zc')
    end
end