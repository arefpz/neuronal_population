%% Written by Aref Pariz
% Find the histogram of spike timing differences for taum1=10ms and "taum"
% at stimulation frequency "omega"

clc;clear;
sec=100;
dt=0.1;
omega=12;
taum=6;
ensemble=10;
folder='data/';
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
IND12=linspace(-150,150,301);
IND21=IND12;
for ens=1:ensemble
    Val12(ens,:)=hist(DT12(1:Nspike(ens),ens),IND12);
    Val21(ens,:)=hist(DT21(1:Nspike(ens),ens),IND21);
end
%%
figure;hold on;
errorbar((IND12),(mean(Val12,1)),std(Val12,1,1),'b')
errorbar((IND21),(mean(Val21,1)),std(Val21,1,1),'r')
xlim([-100 100])
legend('g_{1\rightarrow2}','g_{2\rightarrow1}')

title(['\tau_2=',num2str(taum),'\omega=',num2str(omega)]);