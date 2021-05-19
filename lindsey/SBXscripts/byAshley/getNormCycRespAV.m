function [avgNormTC,avgNormResp,semNormResp] = getNormCycRespAV(...
    normTCEaCycAllCells,normCondition,respwin)

[ncycles,nav] = size(normTCEaCycAllCells);
avgNormTC = zeros(size(normTCEaCycAllCells{1,normCondition},1),ncycles,nav);
avgNormResp = zeros(ncycles,nav);
semNormResp = zeros(ncycles,nav);
for icyc = 1:ncycles
    for iav = 1:nav
        avgNormTC(:,icyc,iav) = mean(normTCEaCycAllCells{icyc,iav},2);
        avgNormResp(icyc,iav) = mean(mean(...
            normTCEaCycAllCells{icyc,iav}(respwin,:),1),2);
        semNormResp(icyc,iav) = ste(mean(...
            normTCEaCycAllCells{icyc,iav}(respwin,:),1),2);
    end
end

end