function [p n1 n2] = svgttest(fn1,fn2);
%set your path to the folder with the .svg file
%fn1 and fn2 are the names of the .svg files you want to compare
%output is p of the ttest and n for the number of data points

svg1 = loadsvg(fn1,0.5,1);
n1 = size(svg1,2);
yvals{1} = zeros(1,n1);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
for i = 1:n1
    yvals{1}(i) = mean(dataObjs(i).YData);
end

svg2 = loadsvg(fn2,0.5,1);
n2 = size(svg2,2);
yvals{2} = zeros(1,n1);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
for i = 1:n1
    yvals{2}(i) = mean(dataObjs(i).YData);
end

[h p] = ttest2(yvals{1},yvals{2});