
f1 = fspecial('average');
areas = ['PM'; 'LM'; 'AL'];
iArea = 1;
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
image = areas(iArea,:);

sum_base = 'G:\users\lindsey\analysisLG\experiments';
list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);
resp_avg_all = zeros(26,6,nexp);
ind = [];
for iexp = 1:nexp;

    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_protocol = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};

    base = 'G:\users\lindsey\analysisLG\active mice';
    outDir = fullfile(base, mouse, date);

    fn_ttest = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    load(fn_ttest);

    siz = size(Info_ttest_mat);
    Info_ttest_mat_long = reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]);
    ttest_smooth = reshape(filter2(f1,Info_ttest_mat_long), [siz(1) siz(2) siz(3)]);

    alphaB = [0.05 0.01 0.002 0.0001 0.00001 8.3333e-007];

    for ialpha = 1:6;
        ttest_mask_test(:,:,ialpha) = min(ttest_smooth,[],3) < alphaB(ialpha);
        ttest_mask_test(1:5,1:end,ialpha) = 0;
        ttest_mask_test(1:end, 1:5,ialpha) = 0;
        ttest_mask_test(1:end, 251:end,ialpha) = 0;
        ttest_mask_test(235:end,1:end,ialpha) = 0;
    end

    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
    stack = readtiff(fn_stack);
    siz = size(stack);
    stack_dF = zeros(siz(1), siz(2), 26);
    siz = size(stack_dF);
    if dir == 2
        start = 1;
        for iCond = 1:25
            stack_dF(:,:,iCond) = mean(stack(:,:,start:start+1),3);
            start= start+1;
        end
        stack_dF(:,:,end) = stack(:,:,end);
    else
        stack_dF = stack;
    end

    ttest_mask_ind = reshape(ttest_mask_test, [siz(1)*siz(2) 6]);
    stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) siz(3)]);
    resp_avg = zeros(26,6);
    for ialpha = 1:6
        resp_avg(:,ialpha) = mean(stack_dF_long(ttest_mask_ind(:,ialpha), :),1);
    end
    fn_SNR = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_SNR.mat']);
    load(fn_SNR)
    if SNR>5
        ind = [ind iexp];
    end
    resp_avg_all(:,:,iexp) = resp_avg;
end

resp_avg_norm = zeros(26,6,nexp);
for iexp = 1:nexp
    for ialpha = 1:6
        resp_avg_norm(:,ialpha,iexp) = resp_avg_all(:,ialpha,iexp)./max(resp_avg_all(:,ialpha,iexp),[],1);
    end
end
resp_avg_all_avg = nanmean(resp_avg_norm(:,:,ind),3);
resp_avg_all_avg_sq = reshape(resp_avg_all_avg(1:25,:),[5 5 6]);
figure; 
for ialpha = 1:6
    subplot(3,2,ialpha)
    imagesq(resp_avg_all_avg_sq(:,:,ialpha));
    title(num2str(alphaB(:,ialpha)))
    colormap(gray);
    colorbar
end
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120118';
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_PM_SNRabove5_alpha_range.pdf']);
print(gcf, '-dpdf', fn_out);
