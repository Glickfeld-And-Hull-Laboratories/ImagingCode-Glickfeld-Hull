%used in deconvolution_Jin. need raw fluorescence and outputs from the deconvolution function.
%generate a GUI with subplots:
%1. Raw fluorescnece with Spike events have red dots over them.
%2. deconvolved signal
function [fig] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,sessions, threshold)
num_components = size(TCave,2);
fig = figure('Visible','off');
set(gcf,'Position',2*[300 300 960 480]);
set(gcf,'PaperPosition',2*[300 300 960 480]);

sld = uicontrol('Style','slider',...
    'Min',1,'Max',num_components,'Value',1,'SliderStep',[1/(num_components-1) 1],...
    'Position',[150 20 800 20],...
    'Callback',@surfzlim); %num_components is number of total images

txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

fig.Visible = 'on';
plot_component(1)

    function surfzlim(source, callbackdata)
        k = source.Value;
        plot_component(round(k))
    end

    function plot_component(k)                                              % each cell
        cla
        
        plotred = TCave(:,k).*spk_logic(:,k);
        plotred(plotred==0) = NaN;
        subplot(2,1,1)
        plot(TCave(:,k)); hold on;
        plot(plotred,'ro','MarkerSize',8); ylabel('TCave');
        xlim([1,5000]);
        lims = ylim;
        text(10,lims(2)-20,['cell' num2str(k)]);
        drawnow; hold off

        subplot(2,1,2)
        plot(kernel(:,k)); hold on;
        xlim([1,5000]);ylabel('deconvolved signal');xlabel('frame');
        drawnow; hold off
        
    end

supertitle(['TC ave and deconvolution ' sessions ' threshold ' num2str(threshold)]); hold off


end


