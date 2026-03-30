function [lf1s] = getlin_sum(vmdata)


%get the direction combinations for every given orientation
for plstim = 1:length(vmdata(2,1).oris)
    ddind = floor(vmdata(2,1).oris(plstim))/30;
    dircomb(1:2,plstim) = [ddind-2;ddind+2];
end
dircomb(find(dircomb<0)) = dircomb(find(dircomb<0))+12;
dircomb(find(dircomb>11)) = dircomb(find(dircomb>11))-12;

dircomb(:,1:length(vmdata(2,1).oris)) = dircomb(:,1:length(vmdata(2,1).oris))+1;



%now lets get the linear combination of the two gratings for each plaid
for j = 1:length(vmdata)
    for plstim = 1:length(vmdata(2,1).oris)
        % First subtract off the baseline.
        ph = rem(vmdata(2,j).oris(plstim),1);
        rotamt = -ph*8000;  % 2 hz TF at 16000 sampling freqeuncy
        rr = circshift(vmdata(1,j).mcyc(dircomb(1,plstim),:)',rotamt)';
        lsim(plstim,:) = rr + vmdata(1,j).mcyc(dircomb(2,plstim),:) - 2*vmdata(1,j).blank;
        ff = fft(lsim(plstim,:));
        %cep_vm(j).ldc(plstim) = ff(1)/length(ff);
        lf1s(j,plstim) = abs(2*ff(2)/length(ff));
        %cep_vm(j).lf1ang(plstim) = angle(ff(2));
    end
end









% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)

