
%Made by anthony so sorry for any bugs! 

RecordingLabel = "2021_03_23_1615_g0";
FilePathnoRec = "D:\Anthony\CellSorting\2021_03_23_1615_g0\";
FilePath = FilePathnoRec + RecordingLabel;
Channel = "8";

JuiceTimesAdjtprimestring = "TPrime -syncperiod=1.000000 -tostream=" + FilePath + "_tcat.imec0.ap.SY_384_6_500.txt -fromstream=1," + FilePath + "_tcat.nidq.XD_" + Channel + "_0_0.txt -events=1," + FilePath + "_tcat.nidq.XD_" + Channel + "_1_0.txt," + FilePathnoRec + "JuiceTimesAdj.txt" ;
AllLicksAdjtprimestring =  "TPrime -syncperiod=1.000000 -tostream=" + FilePath + "_tcat.imec0.ap.SY_384_6_500.txt -fromstream=1," + FilePath + "_tcat.nidq.XD_" + Channel + "_0_0.txt -events=1," + FilePathnoRec + "AllLicks.txt," + FilePathnoRec + "AllLicksAdj.txt" ;
