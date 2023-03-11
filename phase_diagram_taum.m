%% Written by Aref Pariz
% Uses the data in "folder" and find the phase of spikes
clc;clear;
sec=50;
dt=0.1;
range=1000*sec/dt;
sec1=1 * 10000;
t=linspace(0,sec,range);
st=1;en=range;
TAUM=[6 10 14];
omega_s=12;
Color1=[linspace(0.8,1,numel(TAUM));linspace(0,0,numel(TAUM));linspace(0,1,numel(TAUM))];
Color2=[linspace(0,0,numel(TAUM));linspace(0,1,numel(TAUM));linspace(0.8,1,numel(TAUM))];
ensemble=10;
folder='data/';
ntau=0;
figure;hold on;
nr=1;nc=2;nf=0;
IND=linspace(0,2*pi,101);
xs=2*pi*omega_s*t;
phi_s=mod(xs,2*pi);
for taum=TAUM
    ntau=ntau+1;
    P2=zeros(ensemble,numel(IND));P1=P2;
    for ens=1:ensemble
        fname=[folder,'result_omega',num2str(omega_s),'_taum',num2str(taum),'_ens',num2str(ens),'.mat'];
        load(fname,'rho');
        x1=phi_s(rho(1,st:en));
        x2=phi_s(rho(2,st:en));
        P1(ens,:)=hist(x1,IND);
        P2(ens,:)=hist(x2,IND);
    end
    
    %%
    %     nf=nf+1;
    subplot(nr,nc,1);hold on
    plot(IND,mean(P1,1),'Color',Color1(ntau,:))
    plot(IND,mean(P2,1),'Color',Color1(ntau,:))
    %     xlim([-50 50])
    
    subplot(nr,nc,2);hold on
    errorbar((IND12),(mean(Val12,1)),std(Val12,1,1),'Color',Color1(ntau,:))
    errorbar((IND21),(mean(Val21,1)),std(Val21,1,1),'Color',Color2(ntau,:))
    xlim([-50 50])
end
title(['\tau_2=',num2str(taum)]);