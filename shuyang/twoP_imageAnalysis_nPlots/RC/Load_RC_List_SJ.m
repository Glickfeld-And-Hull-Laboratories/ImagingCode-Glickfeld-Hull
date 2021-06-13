function pairings = Load_RC_List_SJ()
RC_expt_list_SJ;

already_loaded = []; pairings = cell(2,0);

loadingTotal = true;
while loadingTotal
    
loadingID = true;
while loadingID
    prompt =  '\nEnter an experimental ID, or enter "Q" to quit: '; 
    response = input(prompt,'s'); ID = str2double(response);
    if ID >=1 && ID <=6 && mod(ID,1) == 0
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

loadingAnimal = true;
while loadingAnimal
    prompt =  'Enter an animal number, or enter "D" when done: '; 
    response = input(prompt,'s'); logical = ismember(expt(ID).mouse,response,'rows'); animal = response;
    if any(logical)
        if ismember(pairings{2,end},find(logical))
            disp(['This animal has already been added to experiment (' num2str(ID) ')'])
        else
            pairings{2,end}(end+1) = find(logical);
        end
    elseif response == 'D'
        loadingAnimal = false; 
    elseif response == 'Q'
        return
    else
        fprintf(['\nCould not find Animal ' response ' in the experiment list. Enter "Q" to quit and restart.\n']) 
    end
end

prompt =  'Would you like to add another experiment set? Y/N: '; response = input(prompt,'s');
    switch response 
        case 'Y'; loadingTotal = true;
        case 'N'; loadingTotal = false;
        otherwise; disp('Please enter "Y" for yes and "N" for no.')
    end
end
end







