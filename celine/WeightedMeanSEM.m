
function [Ave, SEM, STD, Med] = WeightedMeanSEM(x, n, DIM)
%function [Ave, SEM, STD, Med] = WeightedMeanSEM(x, n, DIM)
%
%x is a 1D or 2D matrix of observation values
%n is a 1D vector.  The "weight" of each observation = effective number of independent observations for each individual observation
%
%DIM = 1 (average rows together), or =2 (average columns together)
%  note: DIM can be excluded, particularly if x is a 1D vector
%
%Returns:
% Ave = mean
% SEM = std error of mean
% STD = std deviation
% Med = median

bMedianKeepsTraceIntact = 0;

if isempty(x)
    Ave = NaN;
    SEM = NaN;
    STD = NaN;
    Med = NaN;
    return;
end

try
    DIM;
catch
    DIM = find(size(x)==length(n));
    if length(DIM)>1 
        if length(x)==1 
            %warning not necessary if this is a sacalar
        else
            warning(['Could not determine correct DIM in WeightedMeanSEM, using first of DIM=' num2str(DIM)]);
        end
        DIM = DIM(1);
    end
end
if size(x,DIM)~=length(n)
    error('mismatch in length(n) and size(x,DIM)');
end
if all(size(n) == size(x))
    nn=n;
else
    if DIM==1
        nn = repmat(n(:),1,size(x,2));
    elseif DIM==2
        nn = repmat(n(:)',size(x,1),1);
    else
        error('WeightedMeanSEM code cannot yet handle DIM>2');
    end
end

Ave = sum(x.*nn,DIM)/sum(n);
if DIM==1
    AA = repmat(Ave,size(x,1),1);
elseif DIM==2
    AA = repmat(Ave,1,size(x,2));
end
STD = sqrt( sum((x-AA).^2.*nn,DIM)/(sum(n)-1) );
SEM = STD/sqrt(sum(n));
if DIM==1
    if size(x,1)==1
        Med = x;
    else
        if bMedianKeepsTraceIntact
            [dum,I] = sort(MeanNoNaN(x,2)); %sort waveforms accordingly
            x = x(I,:);
            n = n(I);
            Med = interp1(cumsum(n),x,sum(n)/2); %linearly interpolates to find x(cumsum(n)=0.5sum(n))
        else
            %sort each point independently
            [dum,I] = sort(x, 1);
            for w=1:size(x,2)
                xx(:,w) = x(I(:,w),w);
                nnn(:,w) = nn(I(:,w),w);
            end
            
            Med = zeros(1, size(x,2));
            for w=1:size(x,2)
                Med(w) = interp1(cumsum(nnn(:,w)),xx(:,w),sum(nnn(:,w))/2); %linearly interpolates to find x(cumsum(n)=0.5sum(n))
            end
            
        end
    end
elseif DIM==2
    if size(x,2)==1
        Med = x;
    else
        if bMedianKeepsTraceIntact
            %interp1 only works for DIM=1
            [dum,I] = sort(MeanNoNaN(x,1)); %sort accordingly
            x = x(:,I);
            n = n(I);
            Med = interp1(cumsum(n),x',sum(n)/2)'; %linearly interpolates to find x(cumsum(n)=0.5sum(n))
        else
            error('code not done yet');
        end
    end
end
