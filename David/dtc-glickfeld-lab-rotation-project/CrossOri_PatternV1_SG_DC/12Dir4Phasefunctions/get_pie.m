function [p] = get_pie(ZpZcStruct)

 colors = [
    160, 160, 160;   %light grey
    70, 130, 180;    %blue
    0, 100, 0;       %dark green
    0, 0, 0;         %black
    ] / 255;


all_pds = ZpZcStruct.PDSind_byphase;
all_cds = ZpZcStruct.CDSind_byphase;
total = length(ZpZcStruct.Zp);

for i = 1:4
    pds_percent(i) = length(all_pds{i}) / total;
    cds_percent(i) = length(all_cds{i}) / total;
    unclass_percent(i) = 1-(pds_percent(i)+cds_percent(i));
end

data = [mean(pds_percent), mean(cds_percent), mean(unclass_percent)];
figure
p = piechart(data,Direction="clockwise");
p.StartAngle = 0;
p.LegendVisible = 1;
p.ColorOrder = colors;
p.FaceAlpha = 0.5;
p.Names = [{'Pattern'}, {'Component'}, {'Unclassified'}];
p.EdgeColor = [0.5 0.5 0.5];
p.LabelStyle = 'percent';

