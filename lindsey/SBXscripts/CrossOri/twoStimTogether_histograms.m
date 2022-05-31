ds = 'TwoStimTogether_ExptList';
eval(ds)
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
fn_out = fullfile(LG_base, 'Analysis\2P');
nexp = length(expt);
uselist = [2 4; 3 4; 2 7; 3 7];
stimlist = {'vert' 'vert'; 'horz' 'vert'; 'vert' 'horz'; 'horz' 'horz'};
twolist = [5; 6; 8; 9];
ncond = size(uselist,1);
resp_peak = cell(ncond,3);
norm_resp = cell(1,ncond);
isnorm = [];
nmax = [];
for iexp = 1:nexp
    if strcmp(expt(iexp).driver,'SLC')
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        coFolder = expt(iexp).coFolder;
        nrun = length(coFolder);
        run_str = catRunName(cell2mat(coFolder), nrun);
        fprintf([mouse ' ' date '\n'])
        datemouse = [date '_' mouse];
        datemouserun = [date '_' mouse '_' run_str];
        load(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_respData.mat']))
        resp_ind = cell(1,ncond);
        for i = 1:ncond
            resp_ind{i} = find(sum(h_resp(:,uselist(i,:)),2)==2);
        end
        
        for i = 1:ncond
            if length(resp_ind{i})
                for ii = 1:length(resp_ind{i})
                    if ttest2(resp_cell{uselist(i,1)}(resp_ind{i}(ii),:),resp_cell{uselist(i,2)}(resp_ind{i}(ii),:),'tail','both')
                        if ~swtest(resp_cell{uselist(i,1)}(resp_ind{i}(ii),:)) & ~swtest(resp_cell{uselist(i,2)}(resp_ind{i}(ii),:))
                            norm_resp{i} = [norm_resp{i} resp_ind{i}(ii)];
                            figure;
                            subplot(3,1,1)
                            histogram(resp_cell{uselist(i,1)}(resp_ind{i}(ii),:),[-0.3:0.025:4])
                            vline(mean(resp_cell{uselist(i,1)}(resp_ind{i}(ii),:)))
                            resp_peak{i,1} = [resp_peak{i,1} mean(resp_cell{uselist(i,1)}(resp_ind{i}(ii),:))];
                            title(stimlist{i,1})
                            xlim([-0.3 4])
                            subplot(3,1,2)
                            histogram(resp_cell{uselist(i,2)}(resp_ind{i}(ii),:),[-0.3:0.025:4])
                            vline(mean(resp_cell{uselist(i,2)}(resp_ind{i}(ii),:)))
                            resp_peak{i,2} = [resp_peak{i,2} mean(resp_cell{uselist(i,2)}(resp_ind{i}(ii),:))];
                            xlim([-0.3 4])
                            title(stimlist{i,2})
                            subplot(3,1,3)
                            H = histogram(resp_cell{twolist(i)}(resp_ind{i}(ii),:),[-0.3:0.025:4]);
                            xlim([-0.3 4])
                            
                            if ~swtest(resp_cell{twolist(i)}(resp_ind{i}(ii),:))
                                vline(mean(resp_cell{twolist(i)}(resp_ind{i}(ii),:)))
                                resp_peak{i,3} = [resp_peak{i,3} {mean(resp_cell{twolist(i)}(resp_ind{i}(ii),:))}];
                                title([stimlist{i,1} '-' (stimlist{i,2}) ': norm' ])
                                isnorm = [isnorm 1];
                                nmax = [nmax 1];
                            else
                                vline(H.BinEdges(islocalmax(H.Values,'MinSeparation',3,'MinProminence',2)));
                                resp_peak{i,3} = [resp_peak{i,3} {H.BinEdges(islocalmax(H.Values,'MinSeparation',3,'MinProminence',2))}];
                                title([stimlist{i,1} '-' (stimlist{i,2}) ': not norm' ])
                                isnorm = [isnorm 0];
                                nmax = [nmax length(resp_peak{i,3})];
                            end
                            suptitle([mouse ' ' date ' cell ' num2str(resp_ind{i}(ii))])
                        end
                    end
                end
            end
        end
    end
end