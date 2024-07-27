function [psth] = psthfxn(StimIdent,NeuronByTrialType,NeuronNumber)

% 
% NeuronNum=size(NeuronByTrialType(1).Stim,2)
% TrialNum=size(NeuronByTrialType(1).Stim(1).BLResp,2)

% NeuronNumber=5

bin_size=0.01;
% b=zeros(TrialNum,50);
binedgesbl=0:0.01:.49;
binedgessd=0:0.01:.29;
binedgesbn=0:0.01:.04;
binedgesss=0:0.01:.29;

% Trials=1:25
for h=1:length(StimIdent)
for i=1:size(NeuronByTrialType(h).Stim(NeuronNumber).BLResp,2)
    a=NeuronByTrialType(h).Stim(NeuronNumber).BLResp(i).Stim;
    [trialpsthBL, bins]=hist(a,binedgesbl);
    b(i,:)=trialpsthBL;
    c=NeuronByTrialType(h).Stim(NeuronNumber).SDResp(i).Stim;
    [trialpsthSD, bins]=hist(c,binedgessd);
    d(i,:)=trialpsthSD;
    e=NeuronByTrialType(h).Stim(NeuronNumber).BTNResp(i).Stim;
    [trialpsthBN, bins]=hist(e,binedgesbn);
    f(i,:)=trialpsthBN;
     y=NeuronByTrialType(h).Stim(NeuronNumber).SSResp(i).Stim;
    [trialpsthSS, bins]=hist(y,binedgesss);
    z(i,:)=trialpsthSS;
end

psthbl = sum(b,1) ./size(NeuronByTrialType(h).Stim(NeuronNumber).BLResp,2) ./ bin_size;
psth2bl(h,:) = smoothdata(psthbl,'gaussian');
psthSD = sum(d,1) ./size(NeuronByTrialType(h).Stim(NeuronNumber).BLResp,2) ./ bin_size;
psth2SD(h,:) = smoothdata(psthSD,'gaussian');
psthBN = sum(f,1) ./size(NeuronByTrialType(h).Stim(NeuronNumber).BLResp,2) ./ bin_size;
psth2BN(h,:) = smoothdata(psthBN,'gaussian');
psthSS = sum(z,1) ./size(NeuronByTrialType(h).Stim(NeuronNumber).BLResp,2) ./ bin_size;
psth2SS(h,:) = smoothdata(psthSS,'gaussian');
end
%%
binedgesss=binedgesss+.85;
binedgesbn=binedgesbn+0.8;
binedgessd=binedgessd+0.5;
binedges=[binedgesbl binedgessd binedgesbn binedgesss];
psth=[psth2bl psth2SD psth2BN psth2SS];
%%
for i=1:length(StimIdent)
    figure
plot(binedges,psth(i,:), 'color','b','LineWidth',3)
box off
title(['Neuron ',num2str(NeuronNumber),' Response to ',num2str(StimIdent(i)),' CPD'])
xline([.5 .5])
xline([.85 .85])
hold off
ylabel('Firing Rate (Hz)')
xlabel('Time (S)')
% xticklabels(StimIdent)
end
