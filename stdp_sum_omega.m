%% Written by Aref Pariz
% analyze the simulated data in "folder" and find the spike timing
% differences between couled neurons "dT" and the convolution of STDP and dT
clc;clear;
dt=0.1;
OMEGA=0:2:48;
TAUM=[6 10 14];%4:2:20;
Color1=[linspace(0.8,1,numel(OMEGA));linspace(0,0,numel(OMEGA));linspace(0,1,numel(OMEGA))];
Color2=[linspace(0,0,numel(OMEGA));linspace(0,1,numel(OMEGA));linspace(0.8,1,numel(OMEGA))];
ensemble=10;
folder='data/';
figure;hold on;
nr=3;nc=4;nf=8;
ntaum=0;
for taum=TAUM
    sum_all1=zeros(1,numel(OMEGA));
    sum_all2=zeros(1,numel(OMEGA));
    ntaum=ntaum+1;
    nomega=0;
    for omega=OMEGA
        nomega=nomega+1;
        DT12=[];DT21=[];
        Nspike=zeros(1,ensemble);
        for ens=1:ensemble
            fname=[folder,'result_omega',num2str(omega),'_taum',num2str(taum),'_ens',num2str(ens),'.mat'];
            load(fname,'rho');
            ts=[0 0];n1=0;
            for ii=1:size(rho,2)
                fired=rho(:,ii);
                if sum(fired)
                    ts(fired)=ii;
                    n1=n1+1;
                    DT12(n1,ens)=dt*(ts(2)-ts(1));
                    DT21(n1,ens)=dt*(ts(1)-ts(2));
                end
            end
            Nspike(ens)=n1;
        end
        %%
        xind=101;
        yind=fix(xind/2)+1;
        IND=linspace(-150,150,xind);
        sum_stdp_hist1=zeros(1,ensemble);
        sum_stdp_hist2=zeros(1,ensemble);
        Val12=zeros(ensemble,numel(IND));
        Val21=zeros(ensemble,numel(IND));
        for ens=1:ensemble
            Val12(ens,:)=hist(DT12(1:Nspike(ens),ens),IND);
            Val21(ens,:)=hist(DT21(1:Nspike(ens),ens),IND);
            x1=IND(1:yind);
            x2=IND(yind:end);
            f1=-exp(x1/10);
            f2=exp(-x2/10);
            sum_stdp_hist1(ens)=sum(f1.*Val12(ens,1:yind))+sum(f2.*Val12(ens,yind:end));
            sum_stdp_hist2(ens)=sum(f1.*Val21(ens,1:yind))+sum(f2.*Val21(ens,yind:end));
        end
        %%
        sum_all1(nomega)=mean(sum_stdp_hist1);
        sum_all2(nomega)=mean(sum_stdp_hist2);
    end
    subplot(1,2,1);hold on;
    plot(OMEGA,sum_all1);
    subplot(1,2,2);hold on;
    plot(OMEGA,sum_all2);
    LEG{ntaum}=(['\tau_m=',num2str(taum)]);
end
legend(LEG)