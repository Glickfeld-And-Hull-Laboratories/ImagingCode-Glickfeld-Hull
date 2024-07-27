function [mfr1] = scottmakeraster(Spike,n)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
spks1=Spike;
sd=n;
fr1 = 1000*smooth(spks1, n);
mfr1= mean(fr1,1);
sdfr1=std(fr1,[],1);
sefr1=sdfr1/sqrt(n);

m=2500;
t=[1:m];
e1=zeros(1,m);
e2=zeros(1,m);
e1(1,:)=mfr1(:,1:m)-sefr1(:,1:m);
e2(1,:)=mfr1(:,1:m)+sefr1(:,1:m);

% h2=subplot(2,1,2);
% hold on
% L(1)=fill([t t(end:-1:1)],[e1,e2(end:-1:1)],'g','EdgeColor','g','FaceColor','g');
% 
% plot(t,mfr1,'Color','k','LineStyle','-')
% xlabel('Time (ms)')
% ylabel('firing rate (spk/s)')
%  hold off 
end

