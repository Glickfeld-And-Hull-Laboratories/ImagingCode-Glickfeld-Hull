
function [SOIDV, SOI, BySec, ByMouse] = HistoVTASummaryDV(OUTPUT)
%function [SOIDV, SOI, BySec, ByMouse] = HistoVTASummaryDV(OUTPUT)
%
%OUTPUT = result of HistoOfflineAnalysis_SB3
%H1 and H2 first column = 'virus' from HistoOfflineAnalysis_SB3 (ch 2)
%H1 and H2 second column = 'dye' from HistoOfflineAnalysis_SB3 (ch 1)


load('HistoIX83_Analysis_SB03.mat');
%averaged over whole mouse
ByMouse.BKG = BkgFromHist((sum(cat(3,OUTPUT.Bkg(:).count2D),3)), OUTPUT.binIntensity);
ByMouse.H1  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
ByMouse.H2  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);

NUMROW = size(OUTPUT.Hemisphere(1).Block,1); %number of dorsal-->ventral subdividions

%by section
NSec = length(OUTPUT.Bkg);
BySec.H1  = zeros(NSec, 2);
BySec.H2  = zeros(NSec, 2);

BySecDV.H1  = zeros(NUMROW, 2, NSec);
BySecDV.H2  = zeros(NUMROW, 2, NSec);
for k=1:NSec
    BySec.BKG = ByMouse.BKG; %copy for completeness
    %row, column, section
    BySec.H1(k,:)  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    BySec.H2(k,:)  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    for nR=1:NUMROW
        %row, column, section
        BySecDV.H1(nR,:,k)  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(nR,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
        BySecDV.H2(nR,:,k)  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(nR,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    end
end

%find the strongest 15-sections (using dye)
%6/14/21: CHANGE TO BE MAX VIRUS - more reliable than dye
DyePerSection = BySec.H1(:,2) + BySec.H2(:,2);
WIDTH = 15; %number of sections
for k=1:NSec-WIDTH+1
    TMP(k) = sum(DyePerSection(k:k+WIDTH-1));
end
[~,kBest] = max(TMP);

%Sections Of Interest
SOI.BKG = ByMouse.BKG; %copy for completeness
SOI.SecIndx = kBest:kBest+WIDTH-1;
SOI.H1   = sum(BySec.H1(SOI.SecIndx,:), 1);
SOI.H2   = sum(BySec.H2(SOI.SecIndx,:), 1);

SOIDV.BKG = ByMouse.BKG; %copy for completeness
SOIDV.SecIndx = kBest:kBest+WIDTH-1;
SOIDV.H1   = sum(BySecDV.H1(:,:,SOI.SecIndx), 3);
SOIDV.H2   = sum(BySecDV.H2(:,:,SOI.SecIndx), 3);

set(0,'DefaultFigureWindowStyle','docked')
%figure(1); clf; imagesc(log10(sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,SOI.SecIndx).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('Hemisphere 1');
%figure(2); clf; imagesc(log10(sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,SOI.SecIndx).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('Hemisphere 2');
%figure(3); clf; imagesc(log10(sum(cat(3,OUTPUT.Bkg(:).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('background');
