%% Written by Aref Pariz
% Uses the data in "folder" and take sum over last 10s of synaptic weights
% and plot the average value of synaptic weights for the second neuron's 
% taum and the stimulation frequency
clc
clear
OMEGA=0:2:50;
TAUM=4:1:20; 
ensemble=10;
Gavg12=zeros(numel(OMEGA),numel(TAUM),ensemble);
Gavg21=zeros(numel(OMEGA),numel(TAUM),ensemble);
folder='data/';
n1=0;
for omega=OMEGA
    n1=n1+1;n2=0;
    for taum=TAUM
        n2=n2+1;
        for ens=1:10
            fname=[folder,'result_omega',num2str(omega),'_taum',num2str(taum),'_ens',num2str(ens),'.mat'];
            load(fname,'G12','G21');
            Gavg12(n1,n2,ens+1)=mean(G12(90e4:end));
            Gavg21(n1,n2,ens+1)=mean(G21(90e4:end));
        end
    end
end
%%
figure(1);
subplot(1,2,1)
imagesc(TAUM,OMEGA,mean(Gavg12,3));
subplot(1,2,2)
imagesc(TAUM,OMEGA,mean(Gavg21,3));