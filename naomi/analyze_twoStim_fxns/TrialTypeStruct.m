function [Neuron,Neuron4Raster] = TrialTypeStruct(StimIdent,Trials,NeuronSDR,wholetrial)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
%% Separates all data by Neuron and Trial Type
for x=1:length(StimIdent)
% for x=2
    %% the number of stimulus configurations
c=find(Trials(:,4)==StimIdent(x));
for h=1:size(NeuronSDR,2)
for i=1:length(c)
    %%%% Use Neuron matrix for analysis
Neuron.Neuron(h).Stim(x).BLResp(i)=NeuronSDR(h).BLResp(c(i));
Neuron.Neuron(h).Stim(x).BLHzCount(i)=NeuronSDR(h).BLHzCount(c(i));
Neuron.Neuron(h).Stim(x).BLHz(i)=NeuronSDR(h).BLHz(c(i));
Neuron.Neuron(h).Stim(x).SDResp(i)=NeuronSDR(h).SDResp(c(i));
Neuron.Neuron(h).Stim(x).SDHzCount(i)=NeuronSDR(h).SDHzCount(c(i));
Neuron.Neuron(h).Stim(x).SDHz(i)=NeuronSDR(h).SDHz(c(i));
Neuron.Neuron(h).Stim(x).SSResp(i)=NeuronSDR(h).SSResp(c(i));
Neuron.Neuron(h).Stim(x).SSHzCount(i)=NeuronSDR(h).SSHzCount(c(i));
Neuron.Neuron(h).Stim(x).SSHz(i)=NeuronSDR(h).SSHz(c(i));
%%%% Use Neuron4Raster for rasters
Neuron4Raster.Neuron(h).Stim(x).BLResp(i)=wholetrial(h).BLResp(c(i));
Neuron4Raster.Neuron(h).Stim(x).BLHzCount(i)=wholetrial(h).BLHzCount(c(i));
Neuron4Raster.Neuron(h).Stim(x).BLHz(i)=wholetrial(h).BLHz(c(i));
Neuron4Raster.Neuron(h).Stim(x).SDResp(i)=wholetrial(h).SDResp(c(i));
Neuron4Raster.Neuron(h).Stim(x).SDHzCount(i)=wholetrial(h).SDHzCount(c(i));
Neuron4Raster.Neuron(h).Stim(x).SDHz(i)=wholetrial(h).SDHz(c(i));
end
end
end
end