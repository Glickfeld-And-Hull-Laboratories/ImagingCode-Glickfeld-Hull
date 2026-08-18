% instruction sheet for VIP experiments
instructions.sess_list= [18]; %the sessions to concatenate and analyze
instructions.ds='DART_V1_YM90K_Celine';
instructions.experimentFolder = 'VIP_YM90K';
instructions.refDay='2';
instructions.load_retino = true;
instructions.targetCon = [.125 .25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [20.0000  1000.0000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;