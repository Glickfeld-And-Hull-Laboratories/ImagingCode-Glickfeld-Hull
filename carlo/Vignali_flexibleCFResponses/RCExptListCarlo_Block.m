%{
look at:
CS+ after CS-, vv -- 

if sum(expt{i}.order(1:2) == [1,2])==2
pM = 0; mP = 1;
elseif sum(expt{i}.order(1:2) == [2,1])==2
pM = 1; mP = 0;
end


CS+ before CS-, vv
all CS-
all CS+
%}

expt{1}.name=           'Block'; %Day1
expt{1}.mouse =         {'1081','1081','1081','1081'};
expt{1}.date =          {'2102140','2102141','2102150','2102151'};
expt{1}.run =           {'000','001','001','000'};
expt{1}.ttl =           [1,1,1,1];
expt{1}.noreg =         [0,0,0,0];
expt{1}.imgreglaseron = [0,0,0,0];
expt{1}.whichBlock =    [0,1,0,1];
expt{1}.order =         [1,2,2,1];%1 = first; 2 = second

expt{2}.name=           'Block'; %Day1
expt{2}.mouse =         {'1702','1702','1702','1702'};
expt{2}.date =          {'2103310','2103311','2104010','2104011'};
expt{2}.run =           {'000','000','000','000'};
expt{2}.ttl =           [1,1,1,1];
expt{2}.noreg =         [0,0,0,0];
expt{2}.imgreglaseron = [0,0,0,0];
expt{2}.whichBlock =    [0,1,0,1];
expt{2}.order =         [2,1,1,2];

%{
expt{3}.name=           'Block'; %Day1
expt{3}.mouse =         {'1705','1705','1705','1705'};
expt{3}.date =          {'2105160','2105161','2105170','2105171'};
expt{3}.run =           {'000','000','001','000'};
expt{3}.ttl =           [1,1,1,1];
expt{3}.noreg =         [0,0,0,0];
expt{3}.imgreglaseron = [0,0,0,0];
expt{3}.whichBlock =    [0,1,0,1];
expt{3}.order =         [2,1,1,2];
%}

expt{3}.name=           'Block'; %Day1
expt{3}.mouse =         {'1706','1706','1706','1706'};
expt{3}.date =          {'2106080','2106081','2106090','2106091'};
expt{3}.run =           {'000','001','000','001'};
expt{3}.ttl =           [1,1,1,1];
expt{3}.noreg =         [0,0,0,0];
expt{3}.imgreglaseron = [0,0,0,0];
expt{3}.whichBlock =    [0,1,0,1];
expt{3}.order =         [2,1,1,2];

expt{4}.name=           'Block'; %Day1
expt{4}.mouse =         {'1710','1710','1710','1710'};
expt{4}.date =          {'2107220','2107221','2107230','2107231'};
expt{4}.run =           {'000','001','001','000'};
expt{4}.ttl =           [1,1,1,1];
expt{4}.noreg =         [0,0,0,0];
expt{4}.imgreglaseron = [0,0,0,0];
expt{4}.whichBlock =    [0,1,0,1];
expt{4}.order =         [1,2,2,1];

