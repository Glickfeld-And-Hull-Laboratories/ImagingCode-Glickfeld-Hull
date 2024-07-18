function [pairings] = loadRCList_Carlo(expt,seshType)
%% QUESTIONS ASKED WHILE THIS FUNCTION RUNS + EXPECTED ANSWERS
% 
% Enter an experimental ID, or enter "Q" to quit:
%%% 1-X; X being the size of 'expt' or the number of animals in the experiment list
% Enter the session number (interleaved is 1 for naive/2 for PL, blocked 1-4): 
%%% self-explanatory imo; which value you can select depends on whether you
%%% loaded the _Block.m or _Inter.m list; might combine later idk
% Enter an animal number, or enter "d" when done: 
%%% e.g. 1706
% Enter an animal number, or enter "d" when done: 
%%% almost always type "d" here
% Would you like to add another experiment set? y/n: 
%%% almost always type "n" here
% 
% 
% 
%%
already_loaded = []; pairings = cell(3,0);

loadingTotal = true;
while loadingTotal
    
loadingID = true;
while loadingID
    %
    
    prompt =  '\n Enter an experimental ID, or enter "Q" to quit: '; 
    response = input(prompt,'s'); 
    ID = str2double(response);
    if ID >= 1 && mod(ID,1) == 0 %&& ID <= 6
        if ismember(already_loaded,ID)
            disp(['Experiment ID (' num2str(ID) ') is already loaded. Restart script to include more animals.'])
        else
           loadingID = false; already_loaded(end+1) = ID; pairings{1,end+1} = ID;
        end
    elseif response == 'Q'
        return
    else
        disp('Invalid entry. Please enter the experiement ID as an integer from 1-6.') 
    end
end

if seshType ~= 2
     
 prompt = 'Enter the session number (interleaved is 1 for naive/2 for PL, blocked 1-4): ';
    response = sscanf(input(prompt,'s'),'%d'); 
    logical = (expt{1,ID}.whichBlock(1,response)); 
    blockSesh = response; %expt{1,ID}.whichBlock(1,response);
    if logical==1||logical==0
        if ismember(pairings{3,end},find(logical))
            disp(['This session has already been added to the experiment']);
        else
            pairings{3,end}(end+1) = blockSesh;
        end
    elseif response >= 5
        disp(['Only 4 blocked sessions are performed, if you are looking for interleaved the value is 1']);
    elseif response == 'Q'
        return
    else
        disp(['Value entered does not correspond']);
    end
    
loadingAnimal = true;
while loadingAnimal
    prompt =  'Enter an animal number, or enter "d" when done: '; 
    response = input(prompt,'s'); logical = ismember(expt{1,ID}.mouse{1,1},response,'rows'); animal = response;
    if any(logical)
        if ismember(pairings{2,end},find(logical))
            disp(['This animal has already been added to experiment (' num2str(ID) ')'])
        else
            pairings{2,end}(end+1) = find(logical);
        end
    elseif response == 'd'
        loadingAnimal = false; 
    elseif response == 'Q'
        return
    else
        fprintf(['\nCould not find Animal ' response ' in the experiment list. Enter "Q" to quit and restart.\n']) 
    end
    
end

prompt =  'Would you like to add another experiment set? y/n: '; response = input(prompt,'s');
    switch response 
        case 'y'; loadingTotal = true;
        case 'n'; loadingTotal = false;
        otherwise; disp('Please enter "y" for yes and "n" for no.')
    end


elseif seshType==2
        if sum(expt{ID}.order(1:2) == [1,2])==2
        pM = 0; mP = 1;
        elseif sum(expt{ID}.order(1:2) == [2,1])==2
        pM = 1; mP = 0;
        end
    
 prompt = 'CS+ second [1], CS+ first[2], CS- second [3], CS- first[4]: ';
    response = sscanf(input(prompt,'s'),'%d');
    for block = 1:4
        if expt{1,ID}.whichBlock(1,block) == 1 && expt{1,ID}.order(1,block) == 2 && response==1
            thisBlock = block;
        elseif expt{1,ID}.whichBlock(1,block) == 1 && expt{1,ID}.order(1,block) == 1 && response==2
            thisBlock = block;
        elseif expt{1,ID}.whichBlock(1,block) == 0 && expt{1,ID}.order(1,block) == 2 && response==3
            thisBlock = block;
        elseif expt{1,ID}.whichBlock(1,block) == 0 && expt{1,ID}.order(1,block) == 1 && response==4
            thisBlock = block;
        end
    end
    logical = (expt{1,ID}.whichBlock(1,thisBlock)); blockSesh = thisBlock; %expt{1,ID}.whichBlock(1,response);
    if logical==1||logical==0
        if ismember(pairings{3,end},find(logical))
            disp(['This session has already been added to the experiment']);
        else
            pairings{3,end}(end+1) = blockSesh;
        end
    elseif response >= 5
        disp(['Only 4 blocked sessions are performed, if you are looking for interleaved the value is 1']);
    elseif response == 'Q'
        return
    else
        disp(['Value entered does not correspond']);
    end
    
loadingAnimal = true;
while loadingAnimal
    prompt =  'Enter an animal number, or enter "d" when done: '; 
    response = input(prompt,'s'); logical = ismember(expt{1,ID}.mouse{1,1},response,'rows'); animal = response;
    if any(logical)
        if ismember(pairings{2,end},find(logical))
            disp(['This animal has already been added to experiment (' num2str(ID) ')'])
        else
            pairings{2,end}(end+1) = find(logical);
        end
    elseif response == 'd'
        loadingAnimal = false; 
    elseif response == 'Q'
        return
    else
        fprintf(['\nCould not find Animal ' response ' in the experiment list. Enter "Q" to quit and restart.\n']) 
    end
    
end

prompt =  'Would you like to add another experiment set? y/n: '; response = input(prompt,'s');
    switch response 
        case 'y'; loadingTotal = true;
        case 'n'; loadingTotal = false;
        otherwise; disp('Please enter "y" for yes and "n" for no.')
    end
end
end
end








