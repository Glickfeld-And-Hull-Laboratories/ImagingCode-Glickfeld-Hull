%attempting to find appropriate thresholds
nIC = size(data_tc,1);
cut = 1500;
cut1 = 2500;
cut2 = 3500;
range = [15000 40000];

for ic = 1:nIC
    diff_x = diff(data_tc(:,ic));
    %cut = std(diff_x)*5.85;
    cut = thresh(ic);
    ind = find(diff_x>cut);
    ind(find(diff(ind)==1)+1)=[];   %events must be separated by at least one frame. 
    figure;
    subplot(2,1,1)
    plot(diff_x)
    refline(0,cut)
    title(['cell # ' num2str(ic)]);
    subplot(2,1,2)
    plot(data_tc(:,ic)) 
    hold on
    scatter(ind, repmat(12000,length(ind),1), 'filled', 'r')
    title(['Min = ' num2str(cut)]);
    cut
end 



%%
%messing around to find thresholds
nIC = size(data_tc,2);
data_temp = zeros(1000,floor(size(data_tc,1)/1000),nIC);
for ic = 1:nIC
    start = 1;
    for i = 1:floor(size(data_tc,1)/1000)
        data_temp(:,ic,i) = data_tc(start:start+1000-1,ic);
        start = start+1000;
    end
end


cut = 500;
cut1 = 1500;
cut2 = 2500;
range = [15000 40000];
for ic = 1:nIC
    diff_x = diff(squeeze(data_temp(:,ic,1)),[],1);
    ind = find(diff_x(:,1)>cut);
    ind(find(diff(ind)==1)+1)=[];   %events must be separated by at least one frame. 
    figure;
    subplot(2,3,1)
    for ii = 1:length(ind)
        if and(ind(ii)-5>0, ind(ii)+10<size(data_temp,1))
            plot(data_temp(ind(ii)-5:ind(ii)+10,ic,1),'k')
            hold on
        end
    end
    ylim(range)
    title(['Min = ' num2str(cut)])
    %move to next threshold, cut1
    ind_low = find(diff_x(:,1)<cut1);
    ind_high = find(diff_x(:,1)>cut1);
    ind_close = intersect(ind,ind_low);
    ind_above = intersect(ind,ind_high);
    subplot(2,3,5)
    for ii = 1:length(ind_close)
        if ind_close(ii)-5>0
            if ind_close(ii)+10<1000
                plot(data_temp(ind_close(ii)-5:ind_close(ii)+10,ic,1),'b')
                hold on
            end
        end
    end
    ylim(range)
    subplot(2,3,2)
    for ii = 1:length(ind_above)
        if ind_above(ii)-5>0
            if ind_above(ii)+10<1000
                plot(data_temp(ind_above(ii)-5:ind_above(ii)+10,ic,1),'k')
                hold on
            end
        end
    end
    ylim(range)
    title(['Min = ' num2str(cut1)])
    ind_low = find(diff_x(:,1)<cut2);
    ind_high = find(diff_x(:,1)>cut2);
    ind_close = intersect(ind,ind_low);
    ind_above = intersect(ind,ind_high);
    subplot(2,3,6)
    for ii = 1:length(ind_close)
        if ind_close(ii)-5>0
            if ind_close(ii)+10<1000
                plot(data_temp(ind_close(ii)-5:ind_close(ii)+10,ic,1),'b')
                hold on
            end
        end
    end
    ylim(range)
    subplot(2,3,3)
    for ii = 1:length(ind_above)
        if ind_above(ii)-5>0
            if ind_above(ii)+10<1000
                plot(data_temp(ind_above(ii)-5:ind_above(ii)+10,ic,1),'k')
                hold on
            end
        end
    end
    ylim(range)
    title(['Min = ' num2str(cut2)])
    suptitle(['Cell #' num2str(ic) '; ' num2str(length(ind)) ' events'])
end

%thresholds
thresh = [ones(1,nIC).*cut];
thresh([2 11]) = 1200;
thresh([5 8 14]) = 1000;