function [p comp_use n1 n2 n3] = svganova(fn1,fn2,fn3);
%set your path to the folder with the .svg file
%fn1 and fn2 are the names of the .svg files you want to compare
%output is p of the ttest and n for the number of data points

svg1 = loadsvg(fn1,0.5,1);
n1 = size(svg1,2);
yvals{1} = zeros(n1,1);
group{1} = zeros(n1,1);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
for i = 1:n1
    yvals{1}(i,:) = mean(dataObjs(i).YData);
    group{1}(i,:) = 1;
end

svg2 = loadsvg(fn2,0.5,1);
n2 = size(svg2,2);
yvals{2} = zeros(n2,1);
group{2} = zeros(n2,1);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
for i = 1:n2
    yvals{2}(i,:) = mean(dataObjs(i).YData);
    group{2}(i,:) = 2;
end

svg3 = loadsvg(fn3,0.5,1);
n3 = size(svg3,2);
yvals{3} = zeros(n3,1);
group{3} = zeros(n3,1);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
for i = 1:n3
    yvals{3}(i,:) = mean(dataObjs(i).YData);
    group{3}(i,:) = 3;
end

y_all = [yvals{1};yvals{2};yvals{3}];
group_all = [group{1};group{2};group{3}];
[p tab stats] = anova1(y_all,group_all,'off');

comp = multcompare(stats,'Display','off');
comp_use  = array2table(comp);