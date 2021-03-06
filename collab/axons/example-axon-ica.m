%%
setpath;

global USER_PROFILE
USER_PROFILE = 'lindsey'; % 

protdir = dirs.images
%%
newdir = '2011\March\110317\LG26_OriC';
seqfile 

%% register (single plane)
expt = frGetExpt(newdir);
frRegister(newdir,'Overwrite',true,'Oversampling',10,...
               'DoRecurse',false,'Engine','subpixel');
           
%% compute principal components           
frPrinComp(newdir);

%% load principal components
pcs = load(fullfile(expt.dirs.analrootpn,expt.filenames.pca_usv))
nt = size(pcs.v,1);
tt = [0:nt-1]/frGetFrameRate;

%% show principal components, spatial filters
nCs = 16;
figure;
sm = stackFilter(pcs.U,1.5);
ax=[];
for pc = 1:nCs;
    ax(pc)=subax(4,4,pc);
    imstretch(sm(:,:,pc),[.001 .999],1.5);
   
    text(.8,.1,num2str(pc),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;
colormap gray;

%% show first principal components, spectrogram time course
nfft = 256;
figure;
for pc = 1:nCs;
    ax(pc)=subplot(4,4,pc);
    drawSpectrogram(double(pcs.v(:,pc)),nfft);
    text(.8,.1,num2str(pc),'fontsize',12,'color','k','fontweight','bold','unit','norm');
end

%% show principal components, time courses
figure;
delta = .05;
tcOffsetPlot(tt,pcs.v(:,1:nCs),delta,'k');
xlim([0 50]);
offsetlabel([1:nCs],delta,' ');
xlabel('Time (s)');

%% reconstruct movies from first components
rec = pcs.u(:,1:nCs)*pcs.s(1:nCs,1:nCs)*pcs.v(1:100,1:nCs)';
rec = reshape(rec,[240,256,100]);
stackAnimator(rec);

%% 
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
%     movtm] = CellsortPCA(thisfn, [], 150, 1, fullfile(analdir,'cellsort'),[]);
% 
% figure;[PCuse] = CellsortChoosePCs(thisfn, mixedfilters) ;
%     
% CellsortPlotPCspectrum(thisfn, CovEvals, PCuse) ;

%% compute independent components
PCuse = 1:150;
mu = .5;
nIC = 16;
termtol = 1e-6;
maxrounds = 400;
mixedsig = pcs.v';
mixedfilters = pcs.U;
CovEvals = diag(pcs.s).^2;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);

% f0 = mean(reg,3);
dt = 1/frGetFrameRate;

%% plot independent filters
sel = 1:16;
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
for ic = sel
    subplot(6,6,ind);
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;

%% plot independent signals
figure;
tcOffsetPlot(tt,ica_sig(1:16,:)');
offsetlabel([1:16],delta,' ');
xlim([0 250]);
xlabel('Time (s)');
