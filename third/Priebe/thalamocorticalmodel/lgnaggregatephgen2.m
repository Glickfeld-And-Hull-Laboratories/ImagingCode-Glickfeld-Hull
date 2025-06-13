function [os,osv,ds,Zp,Zc,Zpv,Zcv,ar] = lgnaggregratephgen2(sx,sy,nc,thresh,c50,positions)

% 12 orientations, 120 degree plaids
% Lag offset = 2 deg;
lagoffset = 2;

if ~exist('positions')
    xs = normrnd(3,sx,1,nc);  
    ys = normrnd(3,sy,1,nc);
else
    xs = positions(:,1); 
    ys = positions(:,2);
    nc = size(positions,1);
end

wv = 10;  %grating wavelength

% Lagged drive
xs2 = normrnd(3,sx,1,nc) + lagoffset;
ys2 = normrnd(3,sy,1,nc);

gratingcont = 0.5;


% Single gratings
oris = 0:30:330;
% xs
% ys
% xs2
% ys2

sph = pi/6;  % Temporal offset for lagged input, expressed in phase.  Likely should change to time domain

% We want the amplitude and phase for 1) gratings and 2) plaids for wach
% cell at each orientation

% 8 phase set here for plaid
phh = 0:pi/4:(7*pi/4);

% 4 phase set here for plaid
%phh = 0:pi/2:(3*pi/2);

% sine phases for response
phs = 0:pi/30:(2*pi-pi/30);
rphs = 0:-pi/30:(-2*pi+pi/30);

grph = zeros(length(oris),2*nc);
plph = zeros(length(oris),length(phh),2*nc);

gramp = zeros(length(oris),2*nc);
plamp = zeros(length(oris),length(phh),2*nc);


% Start with gratings:
for oo = 1:length(oris)
    for j=1:nc
        [r,p1] = sin2d([xs(j) ys(j)],wv,0,oris(oo));
        grph(oo,j) = p1;
        gramp(oo,j) = gratingcont;

        % Lagged.  Note that sin2d does not know about direction, so we
        % need to include that by changing phase.
        [r,p1] = sin2d([xs2(j) ys2(j)],wv,0,oris(oo));  % xs2 and ys2 are offset in space
        grph(oo,j+nc) = p1+sph;  % sph here is the offset in time.
        gramp(oo,j+nc) = gratingcont;
    end
end


% Now plaids
for oo = 1:length(oris)
    for j=1:nc
        [r,p1] = sin2d([xs(j) ys(j)],wv,0,oris(oo)+60);
        kk = sin(p1)+i*cos(p1);
        for ph = 1:length(phh)
            [r,p2] = sin2d([xs(j) ys(j)],wv,phh(ph),oris(oo)-60);
            ll = sin(p2)+i*(cos(p2));
            plph(oo,ph,j) = angle(kk+ll);
            plamp(oo,ph,j) = gratingcont*abs(kk+ll);
        end

        % Lagged
        [r,p1] = sin2d([xs2(j) ys2(j)],wv,0,oris(oo)+60);
        p1 = p1+sph;  % phase offset for temporal lag
        kk = sin(p1)+i*cos(p1);
        for ph = 1:length(phh)
            [r,p2] = sin2d([xs2(j) ys2(j)],wv,0+phh(ph),oris(oo)-60);
            p2 = p2+sph; % phase offset for temporal lag
            ll = sin(p2)+i*(cos(p2));
            plph(oo,ph,j+nc) = angle(kk+ll);
            plamp(oo,ph,j+nc) = gratingcont*abs(kk+ll);
        end
    end
end

            
% OK, now we know the phase and amplitude of each LGN cell to the target
% cell.  Now we need to expand to include LGN rectification and then
% cortical threshold, for each orientation

% Half wave rectify with the amplitude (could be contrast) then take mean
% Is mean correct?  Yes, likely because we want to save by nc

plamp = plamp*100;
gramp = gramp*100;

R = 1/.666;

n = 1;

% R = 1/.83333333333;
% c50 = 20;
% n = 1;
% 
R = 1/.909090909090;
c50 = 10;
n = 1;

% We should transform the response amplitudes by contrast:
gramp = R*(gramp.^n)./(c50.^n + gramp.^n);
plamp = R*(plamp.^n)./(c50.^n + plamp.^n);


% Now generate population input response
% To do so we will generate the time series response for each LGN neuron 
% and half-wave rectify, add them up and put in a cortical threshold
for oo = 1:length(oris)
    gresp = zeros(1,length(phs));
    plresp = zeros(length(phh),length(phs));
    if oris(oo)<180   % Question here.  Is this correct for plaids?
        uph = phs;
    else
        uph = rphs;
    end
    for j=1:nc*2
        gresp = gresp + rect(gramp(oo,j)*cos(uph+grph(oo,j))-0); % half-wave rectification
        for k=1:length(phh)
            plresp(k,:) = plresp(k,:) + rect(plamp(oo,k,j)*cos(uph+plph(oo,k,j)));
        end
    end
    gresp = gresp./(nc*2);  

    ff = fft(gresp);
    f1v(oo) = 2*abs(ff(2))/length(ff);
    ff = fft(rect(gresp-thresh));
    f1(oo) = 2*abs(ff(2))/length(ff);


    plresp = plresp./(nc*2); 
    for k=1:length(phh)
        ff = fft(plresp(k,:));
        pf1v(k,oo) =  2*abs(ff(2))/length(ff);
        ff = fft(rect(plresp(k,:)-thresh));
        pf1(k,oo) = 2*abs(ff(2))/length(ff);
    end
end


% For firing rate
for j=1:length(phh)
    [Rp(j), Rc(j)] = plcomptest(f1, pf1(j,:),120);
    Zp(j) = ztrans_right([Rp(j) Rc(j)],12);
    Zc(j) = ztrans_right([Rc(j) Rp(j)],12);
end
% For Vm
for j=1:length(phh)
    [Rpv(j), Rcv(j)] = plcomptest(f1v, pf1v(j,:),120);
    Zpv(j) = ztrans_right([Rpv(j) Rcv(j)],12);
    Zcv(j) = ztrans_right([Rcv(j) Rpv(j)],12);
end


os = composi(oris,f1);
osv = composi(oris,f1v);
ds = compdsi(oris,f1);
dsv = compdsi(oris,f1v);

xs = xs-mean(xs);
ys = ys-mean(ys);
ar = eig(cov([xs;ys]'))';
if length(ar)==1
    ar =[0 ar];
end

