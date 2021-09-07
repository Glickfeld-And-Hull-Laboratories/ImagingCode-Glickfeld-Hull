function [f1,f2,f1ang,projectedf1,f1mat,f2mat] = compcontrastrevf1f2(cyc)
%
% function [f1,f2,f1ang,projectedf1,f1mat,f2mat] = compcontrastrevf1f2(cyc)
%
% Computes F1 and F2 amplitudes.  The key is that a simple cell
% will have an F1 which has a modulation that systamtically shifts
% in timing and amplitude


% First get F1s and F2s for each stimulus phase 
% We know there are 12 phases, starting at index 2
for j=1:12
    ff = fft(cyc(j+1,:));
    f1s(j) = 2*ff(2)/length(ff);
    f2s(j) = 2*ff(3)/length(ff);
end

% Must get angle to project F1 onto
% Instead of fitting, let's just go through all angles (hack):
f1mat = [real(f1s)' imag(f1s)'];
angs = 0:179;
for a = angs;
    ra = a*pi/180;
    % Define rotation matrix
    rm = [cos(ra) -sin(ra); sin(ra) cos(ra)];
    rotatedvals = f1mat*rm;
    ff = fft(rotatedvals(:,1));
    f1amp(a+1) = abs(2*ff(2)/length(ff));
end
[Y,I] = max(f1amp);
f1 = Y;
f1ang = angs(I);
ra = f1ang*pi/180;
rm = [cos(ra) -sin(ra); sin(ra) cos(ra)];
rotatedvals = f1mat*rm;
projectedf1 = rotatedvals(:,1);

% Easy to compute F2
f2 = abs(mean(f2s));
f2mat = [real(f2s)' imag(f2s)'];
