%used in twophoton_imageAnalysis

function plotTC_Jin(tc_avg, sm, reg_mean, ICuse, frGetFrameRate)
%nmask = max(max(mask_final));
% RGBmat = NaN(nmask,3);
%
% c = [.1, .5, .9]; c2 = NaN(1, nmask);
% for cc = 1:3:nmask
%     c2(1,[cc:(cc+2)]) = c;
% end
%
%
%
% for a = 1:size(RGBmat,1);
%     one = a/nmask;
%     two = c2(1,a);
%     three = (nmask-a)/nmask;
%     RGBmat(a,:) = [one, two, three];
% end
% timeCourses = tc_avg;
% delta = double(5*nanmean(nanstd(timeCourses)));
% [nsamples, ncells] = size(timeCourses);
% offsets = repmat(0:delta:(ncells-1)*delta,nsamples,1);

% tt = [0:size(tc_avg,1)-1]/frGetFrameRate;
% figure;
% for ab = 1:size(RGBmat,1);
%     hold on
%     fig = plot(tt,timeCourses(:,ab)+offsets(:,ab),'color',RGBmat(ab,:));
% end
colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0];


    shift  = 0;
    fig = figure;
    for i = ICuse
        plot(tc_avg(:,i)+shift,'Color',colord(mod(i-1,size(colord,1))+1,:)); hold on;
        %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
        %         hline(shift);
        %     end
        shift = shift+4500;  %10000 for 2P
    end
    set(gca,'YTick',4500:4500:shift);
    set(gca,'YTicklabel',(ICuse));
   
    
    
end