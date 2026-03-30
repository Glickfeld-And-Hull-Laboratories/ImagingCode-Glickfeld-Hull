% Input argument resp_ind is the index of cells that you want plotted

function plotZpZc4PhasePopulation(ZpZcStruct, resp_ind, sz, example_cell, plot_type)

    Zp      = ZpZcStruct.Zp;
    Zc      = ZpZcStruct.Zc;
    PDSind  = ZpZcStruct.PDSind_byphase;
    CDSind  = ZpZcStruct.CDSind_byphase;
    
    phases  = [0, 90, 180, 270];

    colors = [
    160, 160, 160;   %light grey
    70, 130, 180;    %blue
    25, 100, 25;       %dark green
    0, 0, 0;         %black
    ] / 255;

    if nargin < 4
        plot_type = 'default';
    end

    if strcmp(plot_type, 'default')
        figure; 
        movegui('center')
        for phaseIdx = 1:4
            for i = 1:4
                subplot(4,4,i + (phaseIdx-1)*4) % Compute correct subplot index
                scatter(Zc(i,resp_ind), Zp(i,resp_ind),sz, colors(1,:),"filled",'MarkerEdgeColor',[1 1 1],'linewidth', 1.5)
                hold on
                scatter(Zc(i,PDSind{phaseIdx}), Zp(i,PDSind{phaseIdx}),sz, colors(2,:),"filled",'MarkerEdgeColor',[1 1 1],'linewidth', 1.5);
                scatter(Zc(i,CDSind{phaseIdx}), Zp(i,CDSind{phaseIdx}),sz, colors(3,:),"filled",'MarkerEdgeColor',[1 1 1],'linewidth', 1.5);
                xlabel('Zc')
                ylabel('Zp')
                ylim([-4 8])
                xlim([-4 8])
                if i == 1
                    title(['Pattern cells at ', num2str(phases(phaseIdx))]);
                end
                plotZcZpBorders
            end
        end
    elseif strcmp(plot_type, 'example')
        figure
        movegui('center')

        new_resp = setdiff(resp_ind, example_cell);
        new_PDS = setdiff(PDSind{1}, example_cell);
        new_CDS = setdiff(CDSind{1}, example_cell);
        for i = 1:4
            figure
            scatter(Zc(i,new_resp), Zp(i,new_resp), sz, colors(1,:),"filled",'MarkerEdgeColor',[1 1 1],'linewidth', 1.5)
            hold on
            scatter(Zc(i,new_PDS), Zp(i,new_PDS), sz, colors(2,:),"filled",'MarkerEdgeColor',[1 1 1],'linewidth', 1.5)
            scatter(Zc(i,new_CDS), Zp(i,new_CDS), sz, colors(3,:),"filled",'MarkerEdgeColor',[1 1 1],'linewidth', 1.5)
            scatter(Zc(i,example_cell), Zp(i,example_cell), sz, colors(2,:),'filled', 's', 'MarkerEdgeColor',[1 1 1],'linewidth', 1.5)
            xlabel('Zc')
            ylabel('Zp')
            ylim([-4 8])
            xlim([-4 8])
            title(['Pattern cells at ', num2str(phases(i))])
            plotZcZpBorders
        end
    end

end

