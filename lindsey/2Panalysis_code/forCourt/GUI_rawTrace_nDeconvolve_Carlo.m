%used in deconvolution_Jin. need raw fluorescence and outputs from the deconvolution function.
%generate a GUI with subplots:
%1. Raw fluorescnece with Spike events have red dots over them.
%2. deconvolved signal

%Edited 220324 by Carlo - added a horizontal line at the baseline of each
%cell as calculated by the value of the first 50 frames of each trial
%averaged into a single value for average baseline activity of each cell
function [fig] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,sessions,threshold,mean_baselineTC,std_baselineTC)
%{

TCave = (nFrames, nDendrites) >> raw TC for each neuron ( dFoF )
kernel = (nFrames, nDendrites) >> tau (rise/decay) component for each TC
spk_logic = (nFrames, nDendrites) >> logical array of spikes (spk_logic in deconvolution and all_events in first derivative)
sessions = 'mouse _ date _ run' string >> for title
threshold = (1,1) >> user-determined threshold from deconvolution/first derivative
singleValue_baselineTC = (1,nDendrites) >> baseline average for each cell 





%}
num_components = size(TCave,2);
fig = figure('Visible','off');
set(gcf,'Position',2*[0 0 960/2 480/2]);
set(gcf,'PaperPosition',2*[0 0 960/2 480/2]);

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
        yline(mean_baselineTC(1,k),'r'); %what should be the average baseline of the first 50 frames across all trials for a cell
        xlim([1,1000]); %roughly 4 trials and the max duration without returning to baseline
        ylim([0 max(TCave(:,k))+500]);
        lims = ylim;
        text(10,lims(2)-20,['cell' num2str(k)]);
        drawnow; hold off

        subplot(2,1,2)
        plot(kernel(:,k)); hold on;
        yline(std_baselineTC(1,k),'m'); 
        yline(std_baselineTC(1,k)*2,'r'); 
        yline(std_baselineTC(1,k)*3,'Color',[255/255 128/255 0/255]);
        yline(std_baselineTC(1,k)*4,'Color',[255/255 199/255 0/255]);
        xlim([1,1000]);ylabel('denoised signal/difference vector');xlabel('frame');
        ylim([0 max(kernel(:,k))+500]);
        drawnow; hold off
        
    end

supertitle(['TC ave and deconvolution ' sessions ' threshold ' num2str(threshold)]); hold off


end


