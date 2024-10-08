clear all
areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3;
    P = 2;
    matrix = 'SF5xTF5';
    image = areas(iArea,:);
    inj = 'V1';
    nON=12;
    nOFF=12;
    Nshuf = 500;
    SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
    TF_vec0 = [1 2 4 8 15];

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    base = 'G:\users\lindsey\analysisLG\active mice';
    anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120116';

    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    mouse_list = [];
    mice = exp_list.mouse_mat;

    all_goodfits = [];

for iexp = 1:nexp;
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_prot = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};

    base = 'G:\users\lindsey\analysisLG\active mice';    
    outDir = fullfile(base, mouse,date);

    fn_lbub = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
    load(fn_lbub);

    for iCell = goodfit
    	all_goodfits = [all_goodfits; 2.^lbub_fits(iCell,5,4) 2.^lbub_fits(iCell,4,4) lbub_fits(iCell,1,4)];
        hold on
    end
end
fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_all_goodfits.mat']);
save(fn_out, 'all_goodfits');


speed = all_goodfits(:,1)./all_goodfits(:,2);
speed_vec0 = TF_vec0./SF_vec0;

speed_avg = zeros(5,2);
dF_avg = zeros(5,2);
n = zeros(5,1);
for ispeed = 1:5
    if ispeed == 1
        ind = find(speed<=speed_vec0(1,ispeed)+0.00001);
    else
        ind = find(speed<=speed_vec0(1,ispeed)+0.00001 & speed>speed_vec0(1,ispeed-1)+0.00001);
    end
    speed_avg(ispeed,1) = mean(speed(ind,:),1);
    speed_avg(ispeed,2) = std(speed(ind,:),1)./sqrt(size(ind,1));
    dF_avg(ispeed,1) = mean(all_goodfits(ind,3),1);
    dF_avg(ispeed,2) = std(all_goodfits(ind,3),1)./sqrt(size(ind,1));
    n(ispeed,1) = length(ind);
end
figure;
ploterr(speed_avg(:,1), dF_avg(:,1), speed_avg(:,2), dF_avg(:,2), 'logx')
hold on
for ispeed = 1:5
    text(speed_avg(ispeed,1), 0.1, num2str(n(ispeed, 1)));
end
ylim([0 .5]);
xlim([0 750]);
xlabel('Speed');
ylabel('dF/F');
end
title('PM Speed vs dF errorbar')
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_PM_Speed_vs_dF_errorbar.pdf']);
print(gcf, '-dpdf', fn_out);

