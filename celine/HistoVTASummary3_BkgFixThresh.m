
function [SOI, BySec, ByMouse] = HistoVTASummary3_BkgFixThresh(OUTPUT)
%function [SOI, BySec, ByMouse] = HistoVTASummary3_BkgFixThresh(OUTPUT)
%
%OUTPUT = result of HistoOfflineAnalysis_SB3
%FIRST OUTPUT = 'virus' from HistoOfflineAnalysis_SB3 (ch 2)
%SECOND OUTPUT = 'dye' from HistoOfflineAnalysis_SB3 (ch 1)

load('HistoIX83_Analysis_SB03.mat');
%averaged over whole mouse
%ByMouse.BKG = BkgFromHist((sum(cat(3,OUTPUT.Bkg(:).count2D),3)), OUTPUT.binIntensity);
[ByMouse.BKG, ~, NumBlobs] = BkgFromHist(sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3), OUTPUT.binIntensity, 1);
if NumBlobs~=1
    warning('BkgFromHist NumBlobs should be 1 for this Fix hack to work');
end
NSec = length(OUTPUT.Bkg);

%Zero-Out any data where Red < 1950;
VirusThresh = 1950
for h=1:2
    for k=1:NSec
        for r=1:size(OUTPUT.Hemisphere(h).Block,1)
            for c=1:size(OUTPUT.Hemisphere(h).Block,2)
                OUTPUT.Hemisphere(h).Block(r,c,k).count2D(OUTPUT.binIntensity<VirusThresh,:) = 0;
            end
        end
    end
end
%BkgFromHist(sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3), OUTPUT.binIntensity, 1);

%ByMouse.BKG = [0 0];
[ByMouse.H1, ByMouse.Pix1]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
[ByMouse.H2, ByMouse.Pix2]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);

%by section
BySec.H1  = zeros(NSec, 2);
BySec.H2  = zeros(NSec, 2);
BySec.Pix1  = zeros(NSec, 1);
BySec.Pix2  = zeros(NSec, 1);
for k=1:NSec
    BySec.BKG = ByMouse.BKG; %copy for completeness
    %row, column, section
    [BySec.H1(k,:), BySec.Pix1(k,:)]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
    [BySec.H2(k,:), BySec.Pix2(k,:)]  = SumFromHist( sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,k).count2D),3), OUTPUT.binIntensity, ByMouse.BKG);
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
SOI.BKG      = ByMouse.BKG; %copy for completeness
SOI.SecIndx = kBest:kBest+WIDTH-1;
SOI.H1   = sum(BySec.H1(SOI.SecIndx,:), 1);
SOI.H2   = sum(BySec.H2(SOI.SecIndx,:), 1);
SOI.Pix1 = sum(BySec.Pix1(SOI.SecIndx,:), 1);
SOI.Pix2 = sum(BySec.Pix2(SOI.SecIndx,:), 1);


set(0,'DefaultFigureWindowStyle','docked')
%figure(1); clf; imagesc(log10(sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,SOI.SecIndx).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('Hemisphere 1');
%figure(2); clf; imagesc(log10(sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,SOI.SecIndx).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('Hemisphere 2');
%figure(3); clf; imagesc(log10(sum(cat(3,OUTPUT.Bkg(:).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('background');
