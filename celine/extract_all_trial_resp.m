function [finalArray] = extract_all_trial_resp(nd,nOri,nCon,nCells,nTrials,resp_matrix)
%extract_all_trial_resp turns the cell array of response window dfof values
%into a single matrix with one row per trial per cell, and designates the
%day and stimulus conditions

finalArray=NaN((nd*nCells*nTrials),5);
counter=1;
for iDay = 1:nd
    for iOri=1:nOri
        for iCon = 1:nCon
            %turn the data for this condition and day into an array
            thisArray1 = cell2mat(resp_matrix(iDay,iOri,iCon));
            thisArray1=max(thisArray1,0);
            %reshape it into long format
            thisArray2=reshape(thisArray1,[1,size(thisArray1,1)*size(thisArray1,2)])';
            %duplicate the data into the second column
            thisArray2(:,2)=thisArray2(:,1);
            %replace the first column with cell ID
            thisArray2(:,1)=repelem((1:nCells),size(thisArray1,1));
            %make a third column for contrast
            thisArray2(:,3)=iCon;
            %make a fourth column for ori
            thisArray2(:,4)=iOri;
            %make a fifth column for day
            thisArray2(:,5)=iDay;
            %add this to the bottom of the existing data
            finalArray(counter:(counter+size(thisArray2,1)-1),:)=thisArray2;
            counter=counter+size(thisArray2,1);
        end
    end
end
end

