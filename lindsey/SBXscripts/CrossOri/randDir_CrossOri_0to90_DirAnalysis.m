avg_resp_dir_match = zeros(nCells,5,2);
avg_resp_dir_match(:,:,1) = avg_resp_dir(:,1:5,1);
temp = circshift(avg_resp_dir(:,:,2),2,2);
avg_resp_dir_match(:,:,2) = temp(:,1:5);

int = 22.5;
component_match = circshift(avg_resp_dir_match(:,:,1),-45./int,2)+circshift(avg_resp_dir_match(:,:,1),45./int,2);
pattern_match = avg_resp_dir_match(:,:,1);
comp_corr_match = zeros(1,nCells);
patt_corr_match = zeros(1,nCells);
comp_patt_corr_match = zeros(1,nCells);
for iCell = 1:nCells
    comp_corr_match(iCell) = triu2vec(corrcoef(avg_resp_dir_match(iCell,:,2),component_match(iCell,:)));
    patt_corr_match(iCell) = triu2vec(corrcoef(avg_resp_dir_match(iCell,:,2),pattern_match(iCell,:)));
    comp_patt_corr_match(iCell) = triu2vec(corrcoef(component_match(iCell,:),pattern_match(iCell,:)));
end
Rp_match = ((patt_corr_match)-(comp_corr_match.*comp_patt_corr_match))./sqrt((1-comp_corr_match.^2).*(1-comp_patt_corr_match.^2));
Rc_match = ((comp_corr_match)-(patt_corr_match.*comp_patt_corr_match))./sqrt((1-patt_corr_match.^2).*(1-comp_patt_corr_match.^2));
Zp_match = (0.5.*log((1+Rp_match)./(1-Rp_match)))./sqrt(1./(5-3));
Zc_match = (0.5.*log((1+Rc_match)./(1-Rc_match)))./sqrt(1./(5-3));

figure;
for iC = 1:36
    subplot(6,6,iC)
    plot([0:22.5:90],avg_resp_dir_match(iC,:,1))
    hold on
    plot([0:22.5:90],avg_resp_dir_match(iC,:,2))
    plot([0:22.5:90],component_match(iC,:))
    title([num2str(iC) '- Zc: ' num2str(chop(Zc_match(iC),2)) '; Zp: ' num2str(chop(Zp_match(iC),2))])
end

figure; 
subplot(2,2,1)
scatter(Zc(resp_ind),Zp(resp_ind))
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
axis square
hold on
plotZcZpBorders
subplot(2,2,2)
scatter(Zc_match(resp_ind),Zp_match(resp_ind))
xlabel('Zc match')
ylabel('Zp match')
ylim([-4 8])
xlim([-4 8])
axis square
hold on
plotZcZpBorders
subplot(2,2,3)
scatter(Zc(resp_ind),Zc_match(resp_ind))
ylim([-4 8])
xlim([-4 8])
axis square
xlabel('Zc')
ylabel('Zc match')
subplot(2,2,4)
scatter(Zp(resp_ind),Zp_match(resp_ind))
ylim([-4 8])
xlim([-4 8])
axis square
xlabel('Zp')
ylabel('Zp match')