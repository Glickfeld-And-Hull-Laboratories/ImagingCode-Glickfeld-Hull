
clear all
close all
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\ppDiffOri_2024';

train=3;

load(['alldataAnalyzedTR',num2str(train)])

not8 = [1:7];
%not8 = [1,7];
%not8 = [1,2,3,5,6,7];
%not8 = [4];
%not8 = [1:7];
%DV = DV.{ opt, opt8, est, JM, PV, PVemp, logit, sum }
NDC=500;
names = fieldnames(DV);
names = names(1:9);
dv=[0:NDC]/NDC;

numdatasets=6;

%usedatasets =[1:numdatasets]; dataset_nm = 'All';
%usedatasets =[3:4]; dataset_nm = '8m';
usedatasets =[1:2,5:6]; dataset_nm = '6s';

pn = fullfile(fn_out, dataset_nm);
if ~exist(pn)
    mkdir(pn)
end

clear AUROC
kk=0;
for dataset=usedatasets 
    kk=kk+1;    
for k=1:length(names)
for j=1:3
    for difficulty=1:5
        idxfp = (DVAll.Y==8 & DVAll.cond==1 & logical(sum(DVAll.dataset==dataset,2)));
        switch difficulty
            case 1                
                not8 = 8;
            case 2
                not8 = [1,7];
            case 3
                not8 = [2,6];
            case 4        
                not8 = [3,5];
            case 5
                not8 = 4;                
        end
        idxcd = (logical(sum(DVAll.Y==not8,2)) & DVAll.cond==j & logical(sum(DVAll.dataset==dataset,2))) ;
        DVtemp = abs(getfield(DVAll,names{k}));
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>=dv);

        AUROC{k}(kk,j,difficulty) = -trapz(FP,CD);
    end
    
end
end
end

xa=[0,22.5,45,67.5,90];
for k=1:length(names)
%    subplot(3,3,k), plot(xa,mean(squeeze(AUROC{k}(:,1,:))),'r*'), hold on
%    plot(xa,mean(squeeze(AUROC{k}(:,2,:))),'bo'), hold off
    figure
    errorbar(xa,mean(squeeze(AUROC{k}(:,1,:))),std(squeeze(AUROC{k}(:,1,:)))/sqrt(length(usedatasets)),'b')
    hold on
    errorbar(xa,mean(squeeze(AUROC{k}(:,2,:))),std(squeeze(AUROC{k}(:,2,:)))/sqrt(length(usedatasets)),'r')
    errorbar(xa,mean(squeeze(AUROC{k}(:,3,:))),std(squeeze(AUROC{k}(:,3,:)))/sqrt(length(usedatasets)),'k')
    hold off
    title([names{k} ' GCaMP' dataset_nm])
    ylabel('AUROC')
    xlabel('Orientation difference')
    axis([0,90,0,1])
    legend('250','750','ctrl')
    print(fullfile(pn, [names{k} 'decoder_GCaMP' dataset_nm '.pdf']),'-dpdf','-bestfit')
end

% 
% 
% figure
% for k=1:length(names)
%     subplot(3,3,k)
%     for j=1:train
%         idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
%         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
%         DVtemp = getfield(DV,names{k});
% %        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
%             DVtemp=abs(DVtemp);
% %        end
%         CD = mean(DVtemp(idxcd)>dv);
%         FP = mean(DVtemp(idxfp)>=dv);
%         scatter(FP,CD), hold on
%         xlabel('FP'), ylabel('CD'), title(names{k})
%     end
% legend('250','750','control')
% end
% 
% 
% figure
% for k=1:length(names);
%     subplot(3,3,k)
%     for j=1:max(train,2)
%         idxfp = (DV.Y==8 & DV.cond==1 & logical(sum(DV.dataset==usedatasets,2)));
%         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
%         DVtemp = getfield(DV,names{k});
% %        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
%             DVtemp=abs(DVtemp);
% %        end
%         CD = mean(DVtemp(idxcd)>dv);
%         FP = mean(DVtemp(idxfp)>=dv);
%         scatter(dv,CD), hold on
%         xlabel('DV'), ylabel('CD'), title(names{k})
%     end
% legend('250','750','control')
% end
% 
% figure
% for k=1:length(names);
%     subplot(3,3,k)
%     for j=1:max(train,2)
%         idxfp = (DV.Y==8 & DV.cond==1 & logical(sum(DV.dataset==usedatasets,2)));
%         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
%         DVtemp = getfield(DV,names{k});
%  %       if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
%             DVtemp=abs(DVtemp);
% %        end
%         CD = mean(DVtemp(idxcd)>dv);
%         FP = mean(DVtemp(idxfp)>=dv);
%         scatter(dv,FP), hold on
%         xlabel('DV'), ylabel('FP'), title(names{k})
%     end
% legend('250','750','control')
% end





% %usedatasets =[2,3,4,7,8,9,11,12];
% usedatasets=[2:4,6:12];
% %usedatasets=1:12;
% clear AUC
% kk=0;
% for l=1:length(usedatasets)
% %    usedatasets=l;
%     kk=kk+1;
% for k=1:length(names)
% for d=1:5
%     switch d
%         case 1
%             not8=8;
%         case 2
%             not8=[1,7];
%         case 3
%             not8=[2,6];
%         case 4 
%             not8=[3,5];
%         case 5
%             not8=4;
%     end
% 
%     for j=1:2
%         idxfp = (DV.Y==8 & DV.cond==1 & DV.dataset==usedatasets(l));
% 
%         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==usedatasets(l)) ;
% 
% %        if(j==1 & d==1)
% %            mean(idxfp==idxcd)
% %        end
%         DVtemp = getfield(DV,names{k});
% %        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
%            DVtemp=abs(DVtemp);
% %        end
%         CD = mean(DVtemp(idxcd)>dv);
%         FP = mean(DVtemp(idxfp)>=dv);
%         AUC(kk,k,d,j)=-trapz(FP,CD);
%         if(j==1 & d==1)
%             AUC(kk,k,d,j)=0.5;
%         end
%     end
% end
% end
% end
% 
% 
% clear p
% for k=1:length(names)
% figure    
%     muAUC=squeeze(mean(AUC(:,k,:,:)));
%     stdAUC=squeeze(std(AUC(:,k,:,:)))/sqrt(size(AUC,1));
% %    plot([0:4]/4*90,muAUC)
%     errorbar([0:4]'/4*90,muAUC(:,1),stdAUC(:,1)), hold on
%     errorbar([0:4]'/4*90,muAUC(:,2),stdAUC(:,2)), hold off
%     xlabel('Test - Adapter (^o)')
%     ylabel('auROC')
%     title(names{k})
% %    axis([0 90 0.4 1])
%     legend('250','750')
%     for d=1:2
%         [h,p(k,d)]=ttest(squeeze(AUC(:,k,d,1)),squeeze(AUC(:,k,d,2)));
%     end
% end
% 
% usedatasets =[2:4,6:12];
% 
% %usedatasets=1:12;
% % clear AUC
% % for l=1:length(usedatasets)
% % %    usedatasets=l;
% % 
% % for k=1:length(names);
% %     subplot(3,3,k)
% %     for j=1:max(train,2)
% %         idxfp = (DV.Y==8 & DV.cond==j & DV.dataset==usedatasets(l));
% %         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==usedatasets(l)) ;
% %         DVtemp = getfield(DV,names{k});
% % %        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
% %             DVtemp=abs(DVtemp);
% % %        end
% %         CD = mean(DVtemp(idxcd)>dv);
% %         FP = mean(DVtemp(idxfp)>dv);
% %         scatter(FP,CD), hold on
% %         xlabel('FP'), ylabel('CD'), title(names{k})
% %     AUC(l,k,j)=-trapz(FP,CD);
% %     end
% % legend('250','750','control')
% % end
% % end
% 
% 
% 
% 
% % 
% % figure
% % for k=1:length(names)
% %     subplot(3,3,k)
% %     hist(squeeze(AUC(:,k,:)))
% %     title(names{k})
% %     mu(k,:)=mean(squeeze(AUC(:,k,:)));
% %     mustd(k,:)=std(squeeze(AUC(:,k,:)))/sqrt(size(AUC,3));
% %     legend('250','750','control')
% %     
% %     for j=1:max(train,2)
% %     for l=1:3
% %         sp=sqrt((mustd(k,j)^2+mustd(k,l)^2)/2);
% %         t(k,j,l)=abs(mu(k,j)-mu(k,l))/sp/(2/size(AUC,3));
% %         p(k,j,l)=1-normcdf(t(k,j,l));
% %     end
% %     end
% % end
% 
% 
% % for l=1:numdatasets
% %     usedatasets=l;
% %     figure
% % for k=1:length(names);
% %     subplot(3,3,k)
% %     for j=1:max(train,2)
% %         idxfp = (DV.Y==8 & DV.cond==1 & DV.dataset==l);
% %         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==l) ;
% %         DVtemp = getfield(DV,names{k});
% %  %       if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
% %             DVtemp=abs(DVtemp);
% % %        end
% %         CD = mean(DVtemp(idxcd)>dv);
% %         FP = mean(DVtemp(idxfp)>dv);
% %         scatter(dv,CD), hold on
% %         xlabel('DV'), ylabel('CD'), title(names{k})
% %     end
% % legend('250','750','control')
% % end
% % end
% % 
% % for l=1:numdatasets
% %     usedatasets=l;
% %     figure
% % for k=1:length(names);
% %     subplot(3,3,k)
% %     for j=1:max(train,2)
% %         idxfp = (DV.Y==8 & DV.cond==1 & DV.dataset==l);
% %         idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==l) ;
% %         DVtemp = getfield(DV,names{k});
% % %        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
% %             DVtemp=abs(DVtemp);
% % %        end
% %         CD = mean(DVtemp(idxcd)>dv);
% %         FP = mean(DVtemp(idxfp)>dv);
% %         scatter(dv,FP), hold on
% %         xlabel('DV'), ylabel('FP'), title(names{k})
% %     end
% % legend('250','750','control')
% % end
% % end
% 
% figure(1)
% figure(2)
% figure(3)
% 
% ntheta=size(DV.pr,2);
% theta=[1:ntheta]/ntheta*180-90;
% figure
% hold on
% symb{1}='o';
% symb{2}='.';
% symb{3}='-';
% dataset=usedatasets;
% for j=1:2
%     idx=(DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==dataset,2)));
% %    prbar=exp(mean(log(DV.pr(idx,:))));
%     prbar=mean(DV.pr(idx,:));
%     prbar=prbar/sum(prbar);
%     plot(theta,prbar,strcat('b',symb{j}))
%     idx=(DV.Y==6 & DV.cond==j & logical(sum(DV.dataset==dataset,2)));
% %    prbar=exp(mean(log(DV.pr(idx,:))));
%     prbar=mean(DV.pr(idx,:));
%     prbar=prbar/sum(prbar);
%     plot(theta,prbar,strcat('r',symb{j}))
%     idx=(DV.Y==2 & DV.cond==j & logical(sum(DV.dataset==dataset,2)));
% %    prbar=exp(mean(log(DV.pr(idx,:))));
%     prbar=mean(DV.pr(idx,:));
%     prbar=prbar/sum(prbar);
%     plot(theta,prbar,strcat('g',symb{j}))
% end
% hold off
% title('Epr:  o = 250, . = 750, - = control')
% 
% ntheta=size(DV.pr,2);
% theta=[1:ntheta]/ntheta*180-90;
% figure
% hold on
% symb{1}='o';
% symb{2}='.';
% symb{3}='-';
% for j=1:2
%     idx=(DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==dataset,2)));
%     prbar=exp(mean(log(DV.pr(idx,:))));
% %    prbar=mean(DV.pr(idx,:));
%     prbar=prbar/sum(prbar);
%     plot(theta,prbar,strcat('b',symb{j}))
%     idx=(DV.Y==6 & DV.cond==j & logical(sum(DV.dataset==dataset,2)));
%     prbar=exp(mean(log(DV.pr(idx,:))));
% %    prbar=mean(DV.pr(idx,:));
%     prbar=prbar/sum(prbar);
%     plot(theta,prbar,strcat('r',symb{j}))
%     idx=(DV.Y==2 & DV.cond==j & logical(sum(DV.dataset==dataset,2)));
%     prbar=exp(mean(log(DV.pr(idx,:))));
% %    prbar=mean(DV.pr(idx,:));
%     prbar=prbar/sum(prbar);
%     plot(theta,prbar,strcat('g',symb{j}))
% end
% hold off
% title('expElogpr: o = 250, . = 750, - = control')
% 
% 
% % compute bais and variance of estimators.
% clear theta bias_est var_est 
% theta(1,1,:)=[0.25,0.5,0.75,-1,-0.75,-0.5,-0.25,0];
% %theta(1,1,:)=[-0.25,-0.5,-0.75,-1,0.75,0.5,0.25,0];
% figure
% gg=gcf;
% gg=gg.Number;
% %usedatasets =[2,3,4,7,8,9,11,12];
% usedatasets =[2,3,4,6,7,8,9,10,11,12];
% DVtemp2=nan(size(DV.Y));
% 
% n=0;
% for kk=usedatasets
%     n=n+1;
%     for j=1:3
%     for k=1:8
%         idx=(DV.Y==k & DV.cond==j & DV.dataset==kk);
% %        esttemp = (getfield(DV,names{3}));
% 
% 
%         bias_est(n,j,k)=angle(mean(exp(sqrt(-1)*DV.est(idx)*pi)))/pi;
%         var_est(n,j,k)=1-mean(cos(pi*(DV.est(idx)-bias_est(n,j,k))));
%         DVtemp2(find(idx))=angle(exp(sqrt(-1)*(DV.est(idx)-bias_est(n,j,k))))/pi;
% 
% 
% %        idxii=DV.Y==k & DV.dataset==kk;
% %        bt=angle(mean(exp(sqrt(-1)*DV.est(idxii)*pi)))/pi;
% %        DVtemp2(find(idx))=angle(exp(sqrt(-1)*(DV.est(idx)-bt)))/pi;
% 
% 
% %        bias_est(n,j,k)=mod(bias_est(n,j,k)-theta(k)+1,2)-1+theta(k);
% %        var_est(n,j,k)=1-mean(cos(pi*(DV.est(idx)-bias_est(n,j,k))));
%         bias_est(n,j,k)=(mean(mod(DV.est(idx)-(theta(k))+1,2))-1+(theta(k)));
% 
%         DVtemp2(find(idx))=mod(DV.est(idx)-bias_est(n,j,k)+1,2)-1;
% %        var_est(n,j,k)=std(mod(DV.est(idx)-(theta(k))+1,2))^2;
%         var_est(n,j,k)=std(mod(DV.est(idx)-bias_est(n,j,k)+1,2))^2;
% %        var_est(n,j,k)=mean((mod(DV.est(idx)-(theta(k))+1,2)-1).^2);
% 
%         bias_JM(n,j,k)=angle(mean(exp(sqrt(-1)*DV.JM(idx)*pi)))/pi;
%         var_JM(n,j,k)=1-mean(cos(pi*(DV.JM(idx)-bias_JM(n,j,k))));
%         bias_JM(n,j,k)=mod(bias_JM(n,j,k)-theta(k)+1,2)-1+theta(k);
% 
%         bias_PV(n,j,k)=angle(mean(exp(sqrt(-1)*DV.PV(idx)*pi)))/pi;
%         var_PV(n,j,k)=1-mean(cos(pi*(DV.PV(idx)-bias_PV(n,j,k))));
%         bias_PV(n,j,k)=mod(bias_PV(n,j,k)-theta(k)+1,2)-1+theta(k);
% 
%         bias_PVemp(n,j,k)=angle(mean(exp(sqrt(-1)*DV.PVemp(idx)*pi)))/pi;
%         var_PVemp(n,j,k)=1-mean(cos(pi*(DV.PVemp(idx)-bias_PVemp(n,j,k))));
%         bias_PVemp(n,j,k)=mod(bias_PVemp(n,j,k)-theta(k)+1,2)-1+theta(k);
% 
% %         temp=mod(DV.JM(idx)-(k/8*2-1)+1,2)-1+(k/8*2-1);        
% %         bias_JM(n,j,k)=mean(temp);
% %         var_JM(n,j,k)=var(temp);
% % 
% %         temp=mod(DV.PV(idx)-(k/8*2-1)+1,2)-1+(k/8*2-1);        
% %         bias_PV(n,j,k)=mean(temp);
% %         var_PV(n,j,k)=var(temp);
% % 
% %         temp=mod(DV.PVemp(idx)-(k/8*2-1)+1,2)-1+(k/8*2-1);        
% %         bias_PVemp(n,j,k)=mean(temp);
% %         var_PVemp(n,j,k)=var(temp);
%     end
%     dbias_est(n,j,:)=(circshift(bias_est(n,j,:)-theta,1)-circshift(bias_est(n,j,:)-theta,-1));
%     dbias_JM(n,j,:)=(circshift(bias_JM(n,j,:)-theta,1)-circshift(bias_JM(n,j,:)-theta,-1))*8;
%     dbias_PV(n,j,:)=(circshift(bias_PV(n,j,:)-theta,1)-circshift(bias_PV(n,j,:)-theta,-1))*8;
%     dbias_PVemp(n,j,:)=(circshift(bias_PVemp(n,j,:)-theta,1)-circshift(bias_PVemp(n,j,:)-theta,-1))*8;
% 
%     FI_est(n,j,:)=(1+dbias_est(n,j,:)).^2./var_est(n,j,:);
%     FI_JM(n,j,:)=(1+dbias_JM(n,j,:)).^2./var_JM(n,j,:);
%     FI_PV(n,j,:)=(1+dbias_PV(n,j,:)).^2./var_PV(n,j,:);
%     FI_PVemp(n,j,:)=(1+dbias_PVemp(n,j,:)).^2./var_PVemp(n,j,:);
% 
%     figure(gg)
%     subplot(2,2,1), scatter(theta,squeeze(bias_est(n,j,:))),title('est')
%     hold on
%     subplot(2,2,2), scatter(theta,squeeze(bias_JM(n,j,:))),title('JM')
%     hold on
%     subplot(2,2,3), scatter(theta,squeeze(bias_PV(n,j,:))),title('PV')
%     hold on
%     subplot(2,2,4), scatter(theta,squeeze(bias_PVemp(n,j,:))),title('PVemp')
%     hold on
% 
% 
%     figure(gg+1)
%     subplot(2,2,1), scatter(theta,squeeze(var_est(n,j,:))),title('est')
%     hold on
%     subplot(2,2,2), scatter(theta,squeeze(var_JM(n,j,:))),title('JM')
%     hold on
%     subplot(2,2,3), scatter(theta,squeeze(var_PV(n,j,:))),title('PV')
%     hold on
%     subplot(2,2,4), scatter(theta,squeeze(var_PVemp(n,j,:))),title('PVemp')
%     hold on
% 
% 
%     figure(gg+2)
%     subplot(2,2,1), scatter(theta,squeeze(FI_est(n,j,:))),title('est')
%     hold on
%     subplot(2,2,2), scatter(theta,squeeze(FI_JM(n,j,:))),title('JM')
%     hold on
%     subplot(2,2,3), scatter(theta,squeeze(FI_PV(n,j,:))),title('PV')
%     hold on
%     subplot(2,2,4), scatter(theta,squeeze(FI_PVemp(n,j,:))),title('PVemp')
%     hold on
% 
% 
%     end
% % end
% end
% 
% 
% 
%     idx1=logical(sum(DV.dataset==usedatasets,2)) & DV.cond == 1;
%     idx2=logical(sum(usedatasets==DV.dataset,2)) & DV.cond == 2;
%     idx3=logical(sum(usedatasets==DV.dataset,2)) & DV.cond == 3;
% 
% 
%    for n=usedatasets
%        idx = DV.dataset==n & DV.cond==3 ;
%        temp=std(DVtemp2(idx));
%        mutemp=mean(DVtemp2(idx));
%        idx=DV.dataset==n ;
%        DVtemp2(find(idx))=(DVtemp2(idx))/temp;
%    end
% 
%    idx1s=idx1& logical(sum(DV.Y==[1,7,2,6],2));
%    idx2s=idx2& logical(sum(DV.Y==[1,7,2,6],2));
%    idx3s=idx3& logical(sum(DV.Y==[1,7,2,6],2));
%    idx1s=idx1& logical(sum(DV.Y==[1,7],2));
%    idx2s=idx2& logical(sum(DV.Y==[1,7],2));
%    idx3s=idx3& logical(sum(DV.Y==[1,7],2));
% 
%    [h,p]=vartest2(DVtemp2(idx1),DVtemp2(idx2))
%    [h,p]=vartest2(DVtemp2(idx2),DVtemp2(idx3))
%    [h,p]=vartest2(DVtemp2(idx1),DVtemp2(idx3))
% 
% 
% 
%     stdcorr(1)=std(DVtemp2(idx1));
%     stdcorr(2)=std(DVtemp2(idx2));
%     stdcorr(3)=std(DVtemp2(idx3));
%     varcorr=stdcorr.^2;
%     varcorrstd(1)=varcorr(1)*sqrt(2/(length(idx1)-1));
%     varcorrstd(2)=varcorr(2)*sqrt(2/(length(idx2)-1));
%     varcorrstd(3)=varcorr(3)*sqrt(2/(length(idx3)-1));
%     stdcorrstd=sqrt(varcorrstd);
%     figure 
%     errorbar(varcorr,varcorrstd,'o')
%     hold on
%     bar(varcorr)
%     hold off
% 
% 
%     for j=1:8
%         stdcorr2(j,1)=std(DVtemp2(idx1&DV.Y==j));
%         stdcorr2(j,2)=std(DVtemp2(idx2&DV.Y==j));
%         stdcorr2(j,3)=std(DVtemp2(idx3&DV.Y==j));
%         for n=usedatasets
%             stdcorr3(n,j,1)=std(DVtemp2(idx1&DV.Y==j&DV.dataset==n));
%             stdcorr3(n,j,2)=std(DVtemp2(idx2&DV.Y==j&DV.dataset==n));
%             stdcorr3(n,j,3)=std(DVtemp2(idx3&DV.Y==j&DV.dataset==n));
%         end
%     end
%     stdcorr3=stdcorr3(usedatasets,:,:);
%     stdcorr3=reshape(stdcorr3,length(usedatasets)*8,3);
% 
% %     for n=usedatasets
% %         stdcorr3(n,1)=std(DVtemp2(idx1&DV.dataset==n));
% %         stdcorr3(n,2)=std(DVtemp2(idx2&DV.dataset==n));
% %         stdcorr3(n,3)=std(DVtemp2(idx3&DV.dataset==n));
% %     end
% %     stdcorr3=stdcorr3(usedatasets,:);
% %     stdcorr2=stdcorr2;
% %     
% % stdcorr2
% figure
% errorbar(mean(stdcorr3),std(stdcorr3)/sqrt(size(stdcorr3,1)))
% 
% figure
% theta=squeeze(theta);
% [m,idx]=sort(theta);
% for j=1:3
% mu=squeeze(mean(bias_est(:,j,idx)));
% stdmu=squeeze(std(bias_est(:,j,idx))/sqrt(size(bias_est,1)));
% errorbar(m/2*180,mu/2*180,stdmu/2*180)
% hold on
% end
% refline(1)
% hold off
% 
% mu12=squeeze(var_est(:,1,:)./var_est(:,2,:));
% mu12=mu12(:);
% mu12=mu12(idx);
% [h,p]=ttest(mu12-1)
% 
% mu23=squeeze(var_est(:,2,:)./var_est(:,3,:));
% mu23=mu23(:);
% mu23=mu23(idx);
% [h,p]=ttest(mu23-1)
% 
% mu13=squeeze(var_est(:,1,:)./var_est(:,3,:));
% mu13=mu13(:);
% mu13=mu13(idx);
% [h,p]=ttest(mu13-1)
% 
% figure
% bar([mean(mu13),mean(mu23),mean(mu12)])
% hold on
% errorbar(1,mean(mu13),std(mu13)/length(mu13))
% errorbar(2,mean(mu23),std(mu23)/length(mu23))
% errorbar(3,mean(mu12),std(mu12)/length(mu12))
% 
% stop
% 
% clear bias_est var_est std_est std_bias_est
% usedatasets =[2,3,4,6,7,8,9,10,11,12];
% 
% for n=[3,4,5]
% %for nn=usedatasets    
%     for j=1:3
%     for k=1:8
% %         switch k
% %             case 1
% %                 thetas=8;
% %             case 2
% %                 thetas=[1,7];
% %             case 3
% %                 thetas=[2,6];
% %             case 4
% %                 thetas=[3,5];
% %             case 5
% %                 thetas = 4;
% %         end
%         thetas=k;
%         idx=(logical(sum(DV.Y'==thetas',1))' & DV.cond==j & logical(sum(DV.dataset'==usedatasets',1))');
% 
%         esttemp = (getfield(DV,names{n}));
%         bias_est(j,k)=(mean(mod(esttemp(idx)*180-(theta(k))*180+180,360))-180+(theta(k))*180)/2;
%         std_bias_est(j,k)=std(mod(esttemp(idx)*180-(theta(k))*180+180,360))/2/sqrt(sum(idx));
%         MSE_est(j,k)=sqrt(mean((mod(esttemp(idx)*180-(theta(k))*180+180,360)-180).^2)/4);
%         std_est(j,k)=std(mod(esttemp(idx)*180-(theta(k))*180+180,360))/2;
% 
% %        bias_est(j,k)= mean(esttemp(idx)*180/2);
% %        std_est(j,k)=std(mod(esttemp(idx)*180-abs(theta(k))*180+180,360))/2/sqrt(sum(idx));
% 
% 
%     end
%     end
% 
%     figure
%     for j=3:-1:1
%     [m,idx]=sort(theta,'ascend');
%     errorbar((theta(idx))*90,bias_est(j,idx),std_bias_est(j,idx))
% %    plot(theta(idx)*90,MSE_est(j,idx));
%     hold on
%     end
%     refline(1)
%     hold off
%     legend('Control','750','250')
%     xlabel('Presented Orientation (^o)')
%     ylabel('Estimated Orientation (^o)')
% %    ylabel('sqrt Mean square error (^o)')
%     title(names(n))
% 
%     figure
% 
%     for j=3:-1:1
%     [m,idx]=sort(theta,'ascend');
%     plot((theta(idx))*90,std_est(j,idx)), hold on
%     end
% 
%     xlabel('Presented Orientation (^o)')
%     ylabel('Standard Deviation of Estimator')
%     legend('Control','750','250')
%     title(names(n))
%     hold off
% 
% end
% 
% 
% figure(gg)
% usedatasets=[1:10];
% 
% subplot(2,2,1),refline(1)
% subplot(2,2,2),refline(1)
% subplot(2,2,3),refline(1)
% subplot(2,2,4),refline(1)
% 
% figure
%     subplot(2,2,1), scatter(theta,squeeze(mean(bias_est(usedatasets,1,:)))-squeeze(theta)),title('est'),hold on
%     subplot(2,2,1), scatter(theta,squeeze(mean(bias_est(usedatasets,2,:)))-squeeze(theta)),title('est')
%     subplot(2,2,1), scatter(theta,squeeze(mean(bias_est(usedatasets,3,:)))-squeeze(theta)),title('est')
% %    refline(1)
%     subplot(2,2,2), scatter(theta,squeeze(mean(bias_JM(usedatasets,1,:)))-squeeze(theta)),title('JM'),hold on
%     subplot(2,2,2), scatter(theta,squeeze(mean(bias_JM(usedatasets,2,:)))-squeeze(theta)),title('JM')
%     subplot(2,2,2), scatter(theta,squeeze(mean(bias_JM(usedatasets,3,:)))-squeeze(theta)),title('JM')
% %    refline(1)
%     subplot(2,2,3), scatter(theta,squeeze(mean(bias_PV(usedatasets,1,:)))-squeeze(theta)),title('PV'),hold on
%     subplot(2,2,3), scatter(theta,squeeze(mean(bias_PV(usedatasets,2,:)))-squeeze(theta)),title('PV')
%     subplot(2,2,3), scatter(theta,squeeze(mean(bias_PV(usedatasets,3,:)))-squeeze(theta)),title('PV')
% %    refline(1)
%     subplot(2,2,4), scatter(theta,squeeze(mean(bias_PVemp(usedatasets,1,:)))-squeeze(theta)),title('PVemp'),hold on
%     subplot(2,2,4), scatter(theta,squeeze(mean(bias_PVemp(usedatasets,2,:)))-squeeze(theta)),title('PVemp')
%     subplot(2,2,4), scatter(theta,squeeze(mean(bias_PVemp(usedatasets,3,:)))-squeeze(theta)),title('PVemp')
% %    refline(1)
% legend('250','750','control')
% 
% figure
%     subplot(2,2,1), scatter(theta,squeeze(mean(FI_est(usedatasets,1,:)))),title('est'),hold on
%     subplot(2,2,1), scatter(theta,squeeze(mean(FI_est(usedatasets,2,:)))),title('est')
%     subplot(2,2,1), scatter(theta,squeeze(mean(FI_est(usedatasets,3,:)))),title('est')
% %    refline(1)
%     subplot(2,2,2), scatter(theta,squeeze(mean(FI_JM(usedatasets,1,:)))),title('JM'),hold on
%     subplot(2,2,2), scatter(theta,squeeze(mean(FI_JM(usedatasets,2,:)))),title('JM')
%     subplot(2,2,2), scatter(theta,squeeze(mean(FI_JM(usedatasets,3,:)))),title('JM')
% %    refline(1)
%     subplot(2,2,3), scatter(theta,squeeze(mean(FI_PV(usedatasets,1,:)))),title('PV'),hold on
%     subplot(2,2,3), scatter(theta,squeeze(mean(FI_PV(usedatasets,2,:)))),title('PV')
%     subplot(2,2,3), scatter(theta,squeeze(mean(FI_PV(usedatasets,3,:)))),title('PV')
% %    refline(1)
%     subplot(2,2,4), scatter(theta,squeeze(mean(FI_PVemp(usedatasets,1,:)))),title('PVemp'),hold on
%     subplot(2,2,4), scatter(theta,squeeze(mean(FI_PVemp(usedatasets,2,:)))),title('PVemp')
%     subplot(2,2,4), scatter(theta,squeeze(mean(FI_PVemp(usedatasets,3,:)))),title('PVemp')
% %    refline(1)
% legend('250','750','control')
% 
% % 
% % 
% % NDC=100;
% % l=0;
% % figure
% % % not8 = [1,7];
% % for k=1:3
% % for j=1:3
% %     l=l+1;
% %     subplot(3,3,l) 
% %     for n=1:NDC
% %         dc(n)=(n-1)/(NDC-1);
% %         fpDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==8)>dc(n));
% % %        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y~=8)>dc(n));
% % %        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
% %         cdDV{k,j}(n) = mean(DV{k,j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfDV{k,j}(n) = ((1-fpDV{k,j}(n))*ndist + cdDV{k,j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpDV{k,j},cdDV{k,j},'-o')
% %     hold on
% %     
% %     for n=1:NDC
% %         fpDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y==8)>dc(n));
% % %        cdDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y~=8)>dc(n));
% %         cdDV8{k,j}(n) = mean(DV8{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfDV8{k,j}(n) = ((1-fpDV8{k,j}(n))*ndist + cdDV8{k,j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpDV8{k,j},cdDV8{k,j},'-x')
% %     
% %     for n=1:NDC
% %         fpthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y==8)>dc(n));
% % %        cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y~=8)>dc(n)*2*pi/4);
% %         cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfthetaDV{k,j}(n) = ((1-fpthetaDV{k,j}(n))*ndist + cdthetaDV{k,j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpthetaDV{k,j},cdthetaDV{k,j},'-+')
% % 
% %     for n=1:NDC
% %         fpsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% % %        cdsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %         cdsumDV{k,j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfsumDV{k,j}(n) = ((1-fpsumDV{k,j}(n))*ndist + cdsumDV{k,j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpsumDV{k,j},cdsumDV{k,j},'-s')
% %     
% %     for n=1:NDC
% %         fplogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y==8)>dc(n));
% % %        cdlogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y~=8)>dc(n));
% %         cdlogit{k,j}(n) = mean(logitDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perflogit{k,j}(n) = ((1-fplogit{k,j}(n))*ndist + cdlogit{k,j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fplogit{k,j},cdlogit{k,j},'-*')
% %     plot(0:1,0:1,'k-')
% %     xlabel('FP')
% %     ylabel('CD')
% %     title(['Train = ',num2str(k),' Test = ',num2str(j)])
% %     legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
% %     hold off
% % 
% % end
% % end
% % 
% % 
% % 
% % NDC=100;
% % l=0;
% % figure
% % for j=1:3
% %     l=l+1;
% %     subplot(3,1,l) 
% %     for n=1:NDC
% %         dc(n)=(n-1)/(NDC-1);
% %         fpDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==8)>dc(n));
% % %        cdDVAllj}(n) = mean(DVAll{j}(data{j}.Y~=8)>dc(n));
% % %        cdDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
% %         cdDVAll{j}(n) = mean(DVAll{j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfDVAll{j}(n) = ((1-fpDVAll{j}(n))*ndist + cdDVAll{j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpDVAll{j},cdDVAll{j},'-o')
% %     hold on
% %     
% %     for n=1:NDC
% %         fpDV8All{j}(n) = mean(DV8All{j}(data{j}.Y==8)>dc(n));
% % %        cdDV8All{j}(n) = mean(DV8All{j}(data{j}.Y~=8)>dc(n));
% %         cdDV8All{j}(n) = mean(DV8All{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfDV8All{j}(n) = ((1-fpDV8All{j}(n))*ndist + cdDV8All{j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpDV8All{j},cdDV8All{j},'-x')
% %     
% %     for n=1:NDC
% %         fpthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y==8)>dc(n));
% % %        cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y~=8)>dc(n));
% %         cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfthetaDVAll{j}(n) = ((1-fpthetaDVAll{j}(n))*ndist + cdthetaDVAll{j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpthetaDVAll{j},cdthetaDVAll{j},'-+')
% % 
% %     for n=1:NDC
% %         fpsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% % %        cdsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %         cdsumDVAll{j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));
% %         perfsumDVAll{j}(n) = ((1-fpsumDVAll{j}(n))*ndist + cdsumDVAll{j}(n)*ntarget)/(ndist+ntarget);        
% %     end
% %     plot(fpsumDVAll{j},cdsumDVAll{j},'-s')
% %     
% %     for n=1:NDC
% %         fplogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y==8)>dc(n));
% % %        cdlogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y~=8)>dc(n));
% %         cdlogitAll{j}(n) = mean(logitDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
% %         ndist=sum(data{j}.Y==8);
% %         ntarget=sum(sum(data{j}.Y==not8));        
% %         perflogitAll{j}(n) = ((1-fplogitAll{j}(n))*ndist + cdlogitAll{j}(n)*ntarget)/(ndist+ntarget);
% %     end
% %     plot(fplogitAll{j},cdlogitAll{j},'-*')
% %     plot(0:1,0:1,'k-')
% %     xlabel('FP')
% %     ylabel('CD')
% %     title(['Train = All',' Test = ',num2str(j)])
% %     legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
% %     hold off
% % 
% % end
% % 
% % 
% % 
% % 
% % hand=figure;
% % l=hand.Number-1;
% % for k=1:3
% %     figure(l+2*k-1)
% %     clf
% %     figure(l+2*k)
% %     clf
% % for j=1:3
% %     idx = find(fpDV{k,j} ~= cdDV{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
% %     title(['Optimal with train = ',num2str(k)])
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
% %     title(['Optimal with train = ',num2str(k)])
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpDV8{k,j} ~= cdDV8{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
% %     %legend('FP 250','FP 500','FP 750')
% %     title('Optimal8')
% %     figure(l+2*k)
% %     subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
% %     title('Optimal8')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDV{k,j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDV{k,j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
% %     title('sumDV')
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
% %     title('sumDV')
% %     %legend('CD 250','CD 500','CD 750')    
% %     
% %     idx = find(fplogit{k,j} ~= cdlogit{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
% %     if(j==3)
% %         legend('FP 250','FP 750','FP Control')    
% %     end
% %     title('weighted sum')
% %     xlabel('decision criterion')
% %     figure(l+2*k)
% %     subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
% %     title('weighted sum')
% %     xlabel('decision criterion')
% %     if(j==3)
% %         legend('CD 250','CD 750','CD Control')    
% %     end
% % end
% % end
% % 
% % 
% % 
% % 
% % 
% % hand=figure;
% % l=hand.Number-1;
% % k=1;
% % figure(l+2*k-1)
% % clf
% % figure(l+2*k)
% % clf
% % for j=1:3
% %     idx = find(fpDVAll{j} ~= cdDVAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
% %     title(['Optimal with train = All'])
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
% %     title(['Optimal with train = All'])
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpDV8All{j} ~= cdDV8All{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
% %     %legend('FP 250','FP 500','FP 750')
% %     title('Optimal8')
% %     figure(l+2*k)
% %     subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
% %     title('Optimal8')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDVAll{j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDVAll{j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDVAll{j}(idx)), hold on
% %     title('sumDV')
% % %    legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDVAll{j}(idx)), hold on
% %     title('sumDV')
% % %    legend('CD 250','CD 500','CD 750')    
% %     
% %     idx = find(fplogitAll{j} ~= cdlogitAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
% %     if(j==3)
% %         legend('FP 250','FP 750','FP Control')
% %     end
% %     title('weighted sum')
% %     xlabel('decision criterion')
% %     figure(l+2*k)
% %     subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
% %     title('weighted sum')
% %     if(j==3)
% %         legend('CD 250','CD 750','CD Control')    
% %     end
% %     xlabel('decision criterion')
% % end
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % hand=figure;
% % l=hand.Number-1;
% % for k=1:3
% %     figure(l+2*k-1)
% %     clf
% %     figure(l+2*k)
% %     clf
% % for j=1:3
% %     idx = find(fpDV{k,j} ~= cdDV{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perfDV{k,j};
% %     subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
% %     title(['Optimal with train = ',num2str(k)])
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
% %     title(['Optimal with train = ',num2str(k)])
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpDV8{k,j} ~= cdDV8{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perfDV8{k,j};
% %     subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
% %     %legend('FP 250','FP 500','FP 750')
% %     title('Optimal8')
% %     figure(l+2*k)
% %     subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
% %     title('Optimal8')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perfthetaDV{k,j};
% %     subplot(3,2,3), scatter(dc(idx),fpthetaDV{k,j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,3), scatter(dc(idx),cdthetaDV{k,j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     dc=perfsumDV{k,j};
% %     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
% %     title('sumDV')
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
% %     title('sumDV')
% %     %legend('CD 250','CD 500','CD 750')    
% %     
% %     idx = find(fplogit{k,j} ~= cdlogit{k,j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     dc=perflogit{k,j};
% %     subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
% %     if(j==3)
% %         legend('FP 250','FP 750','FP Control')    
% %     end
% %     title('weighted sum')
% %     xlabel('Performance')
% %     figure(l+2*k)
% %     subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
% %     title('weighted sum')
% %     if(j==3)
% %         legend('CD 250','CD 750','CD Control')    
% %     end
% %     xlabel('Performance')
% % end
% % end
% % 
% % 
% % 
% % 
% % 
% % hand=figure;
% % l=hand.Number-1;
% % k=1;
% % figure(l+2*k-1)
% % clf
% % figure(l+2*k)
% % clf
% % for j=1:3
% %     idx = find(fpDVAll{j} ~= cdDVAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perfDVAll{j};
% %     subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
% %     title(['Optimal with train = All'])
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
% %     title(['Optimal with train = All'])
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpDV8All{j} ~= cdDV8All{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perfDV8All{j};
% %     subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
% %     %legend('FP 250','FP 500','FP 750')
% %     title('Optimal8')
% %     figure(l+2*k)
% %     subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
% %     title('Optimal8')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perfthetaDVAll{j};
% %     subplot(3,2,3), scatter(dc(idx),fpthetaDVAll{j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,3), scatter(dc(idx),cdthetaDVAll{j}(idx)), hold on
% %     title('theta estimated')
% %     %legend('CD 250','CD 500','CD 750')
% % 
% %     idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)    
% %     dc=perfsumDVAll{j};
% %     subplot(3,2,4), scatter(dc(idx),fpsumDVAll{j}(idx)), hold on
% %     title('sumDV')
% % %    legend('FP 250','FP 500','FP 750')
% %     figure(l+2*k)
% %     subplot(3,2,4), scatter(dc(idx),cdsumDVAll{j}(idx)), hold on
% %     title('sumDV')
% % %    legend('CD 250','CD 500','CD 750')    
% %     
% %     idx = find(fplogitAll{j} ~= cdlogitAll{j});
% %     minidx=min(idx)-1;
% %     maxidx=max(idx)+1;
% %     idx=[max(minidx,1):min(maxidx,length(dc))];
% %     figure(l+2*k-1)
% %     dc=perflogitAll{j};
% %     subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
% %     if(j==3)
% %         legend('FP 250','FP 750','FP Control')
% %     end
% %     xlabel('performance')
% %     title('weighted sum')
% %     figure(l+2*k)
% %     subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
% %     title('weighted sum')
% %     if(j==3)
% %         legend('CD 250','CD 750','CD Control')    
% %     end
% %     xlabel('performance')
% % end
% % 
% % 
