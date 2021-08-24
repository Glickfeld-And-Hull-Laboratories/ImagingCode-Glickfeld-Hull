%used in twophoton_imageAnalysis

function sm_dsum = plotMask_Jin(sm, plot_mask)

nmask = size(sm,3);
sm_d = zeros(size(sm,1),size(sm,2),nmask,3);
sm_dsum = zeros(size(sm,1),size(sm,2),3);
c = [.1, .5, .9]; c2 = NaN(1, nmask);
for cc = 1:3:nmask
    c2(1,[cc:(cc+2)]) = c;
end
%tc_avg = zeros(size(img_reg,3),nmask);
for ic=1:nmask
    [i,j] = find(sm(:,:,ic)==1);                             %>mean([max(prctile(sm(:,:,ic),96,1)) max(prctile(sm(:,:,ic),96,2))]));
    %tc = zeros(size(i,1),size(img_reg,3));
    for a = 1:size(i,1)
        sm_d(i(a),j(a),ic,1) = ic/nmask;
        sm_d(i(a),j(a),ic,2) = c2(1,ic);
        sm_d(i(a),j(a),ic,3) = (nmask-ic)/nmask;
        %tc(a,:) = (img_reg(i(a),j(a),:));
    end
    % sorttc = sort(tc,2);                                    %Finds a baseline F for each mask by taking the average of the 200 lowest F values for each pixel in the mask.
    %tcbaseF = mean(mean(tc(1:size(i,1),1:200)));
    %tc_avg(:,ic) = mean(tc,1)';
    %     tc_avg(:,ic) = (tc_avg(:,ic)./tcbaseF)-1;
    for a = 1:3
        sm_dsum(:,:,a) = sm_dsum(:,:,a) + sm_d(:,:,ic,a);
    end
end

% sm_dsum2 = sm_dsum;
% for zz = 1:size(sm_dsum,3);
%     [xx,yy] = find(sm_dsum(:,:,zz)>1);
%     for aa = 1:size(xx,1);
%         sm_dsum2(xx(aa),yy(aa),zz) = 1;
%     end
% end
if plot_mask ==1
figure;
fig = image(sm_dsum);
axis image;
end
%if saveData == 1
 %   saveas(fig, [out_dir, 'mask.fig']);
  %  print([out_dir, 'mask.eps'],'-depsc')
%end
end