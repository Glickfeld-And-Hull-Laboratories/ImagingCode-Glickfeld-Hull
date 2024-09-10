
%%%%%%%%%%%%THE ABOVE SHOULD BE CROSSVALIDATED TO MAKE COMPARISON FAIR>>>>>
clear all
close all

filename{1} = '240529_i1397_runs-002';%_testResp_adapt.mat';
filename{2} = '240529_i1395_runs-002';%_testResp_adapt.mat';  
filename{3} = '240527_i1398_runs-002';%_testResp_adapt.mat';  
filename{4} = '240527_i1396_runs-002';%_testResp_adapt.mat';  
filename{5} = '240523_i1394_runs-002';%_testResp_adapt.mat';  
filename{6} = '240523_i1393_runs-003-004';%_testResp_adapt.mat';  

PCmax=15;
NRmax=PCmax+1;

train=3;

maxiters=2000;
%DV.opt, DV.opt8, DV.est, DV.logit, DV.sum
DV.opt=[];
DV.opt8=[];
DV.est=[];
DV.JM=[];
DV.PV=[];
DV.PVemp=[];
DV.logit=[];
DV.sum=[];
DV.max=[];
DV.dataset=[];
DV.Y=[];
DV.cond=[];
DV.pr=[];
DVAll.opt=[];
DVAll.opt8=[];
DVAll.est=[];
DVAll.JM=[];
DVAll.PV=[];
DVAll.PVemp=[];
DVAll.logit=[];
DVAll.sum=[];
DVAll.max=[];
DVAll.dataset=[];
DVAll.Y=[];
DVAll.cond=[];
DVAll.pr=[];
for n=1:length(filename)

    load([filename{n},'_newFits.mat']);
%    load([filename{n},'_fits'])
    [m,loc]=max(ori_fit);
    
    prefs{n}=loc/180*2*pi;
%    if(exist('theta_90','var'))
        idxn{n}=find(theta_90(3,:)<22.5);% & abs(prefs{n})<pi);
       ['Dataset ',num2str(n),' has ',num2str(length(idxn{n})),' good units using theta_90!']
%    else
%        idxn=find(max(ori_fit)-min(ori_fit)>0.05*max(max(ori_fit)));
%       ['Dataset ',num2str(n),' has ',num2str(length(idxn)),' good units using Amplitude!']
%    end            
    prefs{n}=prefs{n}(idxn{n});
    f{n}=ori_fit(:,idxn{n});
    kappa{n}=abs(fft(log(f{n})));
    kappa{n}=kappa{n}(end,:);
    

    
    for j=1:max(train,2)
        data{j}.X=[];
        data{j}.Y=[];

        for k=1:8   % uses one neurons in ind
            data{j}.X = [data{j}.X;ppResp{j,k}';];
            data{j}.Y = [data{j}.Y;k*ones(size(ppResp{j,k},2),1);];
        end        
        idx=~any(isnan(data{j}.X),2);        
        data{j}.X=data{j}.X(idx,:);
        data{j}.Y=data{j}.Y(idx,1);
        
        
        data{j}.X=data{j}.X(:,idxn{n});
        data{j}.Xraw=[data{j}.X];
    end
%    muX=mean([data{1}.X;data{2}.X;data{3}.X]);
%    stdX=std([data{1}.X;data{2}.X;data{3}.X]);
%    CC=cov([data{1}.X;data{2}.X;data{3}.X]);
    muX=mean([data{1}.X;data{2}.X]);
    stdX=std([data{1}.X;data{2}.X]);
    CC=cov([data{1}.X;data{2}.X]);
%     muX=mean(data{3}.X);
%     stdX=std(data{3}.X);
%     CC=cov(data{3}.X);
    
    [V,D]=eig(CC,'vector');
    
    PCs=min(PCmax,size(data{j}.X,2));
    D=D(max(size(data{1}.X,2)-PCs+1,1):end);
    V=V(:,max(size(data{1}.X,2)-PCs+1,1):end);
    
    
%    D=stdX';
%    V=eye(length(D));
    dataAll.X=[];
    dataAll.Xraw=[];
    dataAll.Y=[];
    dataAll.cond=[];
    for j=1:max(train,2)
        %Zscore or not
%        data{j}.Xraw=bsxfun(@plus,data{j}.Xraw,-muX);
%        data{j}.Xraw=bsxfun(@times,data{j}.Xraw,1./stdX);
%        data{j}.X=bsxfun(@plus,data{j}.X,-muX);
%        data{j}.X=bsxfun(@times,data{j}.X,1./stdX);
        data{j}.X=data{j}.X*V*diag(1./sqrt(D));
        data{j}.X=data{j}.X(:,max(size(data{j}.X,2)-PCs+1,1):end);        
        data{j}.X=[data{j}.X,ones(size(data{j}.X,1),1)];
        
%        if(j<3)
            dataAll.X=[dataAll.X;data{j}.X;];
            dataAll.Xraw=[dataAll.Xraw;data{j}.Xraw;];
            dataAll.Y=[dataAll.Y;data{j}.Y;];
            dataAll.cond=[dataAll.cond;j*ones(size(data{j}.Y));];
%        end
    end

%     % find 20 best units.
%     mAll = VonMesisRegression(size(dataAll.X,2),NRmax);
%     mAll.fit(dataAll.Y/8*2*pi,dataAll.X,maxiters);
%     
%     [m,mind{n}] = sort(mAll.a.mean.^2+mAll.b.mean.^2,'descend'); 
% %    mind{n} = find(abs(mAll.a.mean)./sqrt(diag(mAll.a.var))>2 | abs(mAll.b.mean)./sqrt(diag(mAll.b.var))>2);
%     
%     
%     NR=min(length(mind{n}),NRmax);
%     
%     for j=1:max(train,2)
%         data{j}.X=data{j}.X(:,mind{n}(1:NR));
%     end
% %    clear dataAll mAll
    
    [DVtemp,DVAlltemp,Btemp,BAlltemp]=getDVs(data,dataAll,0,prefs{n}(:,:),kappa{n}(:,:),f{n}(:,:));  % we can use getfield, setfield, and 
                                        % fieldnames to make all DV's are accounted for
                                    
    for j=1:max(train,2)
        Btemp{j}=Btemp{j}(1:end-1,1);
        B{n,j}=V*diag(sqrt(1./D))*Btemp{j};
    end
    
    
% Remove Bias from estimators...

    for j=1:max(train,2)
        DV.opt=[DV.opt;DVtemp{train,j}.opt;];
        DV.opt8=[DV.opt8;DVtemp{train,j}.opt8;];
        DV.est=[DV.est;DVtemp{train,j}.est;];
        DV.JM=[DV.JM;DVtemp{train,j}.JM;];
        DV.PV=[DV.PV;DVtemp{train,j}.PV;];
        DV.PVemp=[DV.PVemp;DVtemp{train,j}.PVemp;];
        DV.logit=[DV.logit;DVtemp{train,j}.logit];
        DV.sum=[DV.sum;DVtemp{train,j}.sum];
        DV.max=[DV.max;DVtemp{train,j}.max];
        DV.dataset=[DV.dataset;n*ones(size(DVtemp{train,j}.opt));];    
        DV.Y=[DV.Y;data{j}.Y;];
        DV.cond=[DV.cond;j*ones(size(data{j}.Y));];
        DV.pr=[DV.pr;DVtemp{train,j}.pr;];
        PVempPrefs = DVtemp{train,train}.PVempPrefs;
        
        DVAll.opt=[DVAll.opt;DVAlltemp{j}.opt;];
        DVAll.opt8=[DVAll.opt8;DVAlltemp{j}.opt8;];
        DVAll.est=[DVAll.est;DVAlltemp{j}.est;];
        DVAll.JM=[DVAll.JM;DVAlltemp{j}.JM;];
        DVAll.PV=[DVAll.PV;DVAlltemp{j}.PV;];
        DVAll.PVemp=[DVAll.PVemp;DVAlltemp{j}.PVemp;];
        DVAll.logit=[DVAll.logit;DVAlltemp{j}.logit];
        DVAll.sum=[DVAll.sum;DVAlltemp{j}.sum];
        DVAll.max=[DVAll.max;DVAlltemp{j}.max];
        DVAll.dataset=[DVAll.dataset;n*ones(size(DVAlltemp{j}.opt));];    
        DVAll.Y=[DVAll.Y;data{j}.Y;];
        DVAll.cond=[DVAll.cond;j*ones(size(data{j}.Y));];
        DVAll.pr=[DVAll.pr;DVAlltemp{j}.pr;];
    end
    
end

%START PLOTTING>>>>>
clear n j k maxiters NR ppResp stdX muX maxiters ans data DVAlltemp DVtemp ind idx CC D loc m ori_fit PCs R_square theta_90 V dataAll
save(['alldataAnalyzedTR',num2str(train)])


% 
% 
% 
% NDC=100;
% l=0;
% figure
% not8 = [1,2,6,7];
% % not8 = [1,7];
% for k=1:3
% for j=1:3
%     l=l+1;
%     subplot(3,3,l) 
%     for n=1:NDC
%         dc(n)=(n-1)/(NDC-1);
%         fpDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==8)>dc(n));
% %        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y~=8)>dc(n));
% %        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
%         cdDV{k,j}(n) = mean(DV{k,j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDV{k,j}(n) = ((1-fpDV{k,j}(n))*ndist + cdDV{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDV{k,j},cdDV{k,j},'-o')
%     hold on
%     
%     for n=1:NDC
%         fpDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y==8)>dc(n));
% %        cdDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y~=8)>dc(n));
%         cdDV8{k,j}(n) = mean(DV8{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDV8{k,j}(n) = ((1-fpDV8{k,j}(n))*ndist + cdDV8{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDV8{k,j},cdDV8{k,j},'-x')
%     
%     for n=1:NDC
%         fpthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y==8)>dc(n));
% %        cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y~=8)>dc(n)*2*pi/4);
%         cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfthetaDV{k,j}(n) = ((1-fpthetaDV{k,j}(n))*ndist + cdthetaDV{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpthetaDV{k,j},cdthetaDV{k,j},'-+')
% 
%     for n=1:NDC
%         fpsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %        cdsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         cdsumDV{k,j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfsumDV{k,j}(n) = ((1-fpsumDV{k,j}(n))*ndist + cdsumDV{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpsumDV{k,j},cdsumDV{k,j},'-s')
%     
%     for n=1:NDC
%         fplogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y==8)>dc(n));
% %        cdlogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y~=8)>dc(n));
%         cdlogit{k,j}(n) = mean(logitDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perflogit{k,j}(n) = ((1-fplogit{k,j}(n))*ndist + cdlogit{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fplogit{k,j},cdlogit{k,j},'-*')
%     plot(0:1,0:1,'k-')
%     xlabel('FP')
%     ylabel('CD')
%     title(['Train = ',num2str(k),' Test = ',num2str(j)])
%     legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
%     hold off
% 
% end
% end
% 
% 
% 
% NDC=100;
% l=0;
% figure
% for j=1:3
%     l=l+1;
%     subplot(3,1,l) 
%     for n=1:NDC
%         dc(n)=(n-1)/(NDC-1);
%         fpDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==8)>dc(n));
% %        cdDVAllj}(n) = mean(DVAll{j}(data{j}.Y~=8)>dc(n));
% %        cdDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
%         cdDVAll{j}(n) = mean(DVAll{j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDVAll{j}(n) = ((1-fpDVAll{j}(n))*ndist + cdDVAll{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDVAll{j},cdDVAll{j},'-o')
%     hold on
%     
%     for n=1:NDC
%         fpDV8All{j}(n) = mean(DV8All{j}(data{j}.Y==8)>dc(n));
% %        cdDV8All{j}(n) = mean(DV8All{j}(data{j}.Y~=8)>dc(n));
%         cdDV8All{j}(n) = mean(DV8All{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDV8All{j}(n) = ((1-fpDV8All{j}(n))*ndist + cdDV8All{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDV8All{j},cdDV8All{j},'-x')
%     
%     for n=1:NDC
%         fpthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y==8)>dc(n));
% %        cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y~=8)>dc(n));
%         cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfthetaDVAll{j}(n) = ((1-fpthetaDVAll{j}(n))*ndist + cdthetaDVAll{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpthetaDVAll{j},cdthetaDVAll{j},'-+')
% 
%     for n=1:NDC
%         fpsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %        cdsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         cdsumDVAll{j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfsumDVAll{j}(n) = ((1-fpsumDVAll{j}(n))*ndist + cdsumDVAll{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpsumDVAll{j},cdsumDVAll{j},'-s')
%     
%     for n=1:NDC
%         fplogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y==8)>dc(n));
% %        cdlogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y~=8)>dc(n));
%         cdlogitAll{j}(n) = mean(logitDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));        
%         perflogitAll{j}(n) = ((1-fplogitAll{j}(n))*ndist + cdlogitAll{j}(n)*ntarget)/(ndist+ntarget);
%     end
%     plot(fplogitAll{j},cdlogitAll{j},'-*')
%     plot(0:1,0:1,'k-')
%     xlabel('FP')
%     ylabel('CD')
%     title(['Train = All',' Test = ',num2str(j)])
%     legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
%     hold off
% 
% end
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% for k=1:3
%     figure(l+2*k-1)
%     clf
%     figure(l+2*k)
%     clf
% for j=1:3
%     idx = find(fpDV{k,j} ~= cdDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8{k,j} ~= cdDV8{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogit{k,j} ~= cdlogit{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')    
%     end
%     title('weighted sum')
%     xlabel('decision criterion')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
%     title('weighted sum')
%     xlabel('decision criterion')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
% end
% end
% 
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% k=1;
% figure(l+2*k-1)
% clf
% figure(l+2*k)
% clf
% for j=1:3
%     idx = find(fpDVAll{j} ~= cdDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8All{j} ~= cdDV8All{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogitAll{j} ~= cdlogitAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')
%     end
%     title('weighted sum')
%     xlabel('decision criterion')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
%     title('weighted sum')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
%     xlabel('decision criterion')
% end
% 
% 
% 
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% for k=1:3
%     figure(l+2*k-1)
%     clf
%     figure(l+2*k)
%     clf
% for j=1:3
%     idx = find(fpDV{k,j} ~= cdDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDV{k,j};
%     subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8{k,j} ~= cdDV8{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDV8{k,j};
%     subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfthetaDV{k,j};
%     subplot(3,2,3), scatter(dc(idx),fpthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx),cdthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     dc=perfsumDV{k,j};
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogit{k,j} ~= cdlogit{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     dc=perflogit{k,j};
%     subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')    
%     end
%     title('weighted sum')
%     xlabel('Performance')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
%     title('weighted sum')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
%     xlabel('Performance')
% end
% end
% 
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% k=1;
% figure(l+2*k-1)
% clf
% figure(l+2*k)
% clf
% for j=1:3
%     idx = find(fpDVAll{j} ~= cdDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDVAll{j};
%     subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8All{j} ~= cdDV8All{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDV8All{j};
%     subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfthetaDVAll{j};
%     subplot(3,2,3), scatter(dc(idx),fpthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx),cdthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     dc=perfsumDVAll{j};
%     subplot(3,2,4), scatter(dc(idx),fpsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx),cdsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogitAll{j} ~= cdlogitAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perflogitAll{j};
%     subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')
%     end
%     xlabel('performance')
%     title('weighted sum')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
%     title('weighted sum')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
%     xlabel('performance')
% end
% 
% 
