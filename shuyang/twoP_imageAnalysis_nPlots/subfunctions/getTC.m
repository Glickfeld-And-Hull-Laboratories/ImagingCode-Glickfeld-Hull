function tc_avg = getTC(img_reg, mask, nmask)

tc_avg = zeros(size(img_reg,3),nmask);
for ic=1:nmask
    %     tc = bsxfun(@times, img_reg, cast(sm(:,:,ic), 'like', img_reg));
    %     tc_avg(:,ic) = mean(mean(tc));
    [i,j] = find(mask(:,:,ic)==1);                             %>mean([max(prctile(sm(:,:,ic),96,1)) max(prctile(sm(:,:,ic),96,2))]));
    tc = zeros(size(i,1),size(img_reg,3));
    for a = 1:size(i,1);
        tc(a,:) = (img_reg(i(a),j(a),:));
    end
    sorttc = sort(tc,2);                                    %Finds a baseline F for each mask by taking the average of the 200 lowest F values for each pixel in the mask.
%     tcbaseF = mean(mean(tc(1:size(i,1),1:200)));
    tc_avg(:,ic) = mean(tc,1)';
    %tc_avg(:,b) = (tc_avg(:,b)./tcbaseF)-1;
end
end