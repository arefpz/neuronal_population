%% 
% The code will plot the distribution of synaptic weights. To use this, you
% should have run the code analyser_taum_less_more.m

figure;
for omega=1:numel(OMEGA)
    subplot(nr,nc,2*nc+1:2*nc+2);hold on;
    plot3(INDEE,OMEGA(omega)*ones(1,numel(INDEE)),mean(GEET{omega}(1:ensemble,:),1),'r');
    plot3(INDEE,OMEGA(omega)*ones(1,numel(INDEE)),mean(GEET{omega}(ensemble+1:2*ensemble,:),1),'b');
    
    subplot(nr,nc,2*nc+3:2*nc+4);hold on;
    plot3(INDEI,OMEGA(omega)*ones(1,numel(INDEI)),mean(GEIT{omega}(1:ensemble,:),1),'r');
    plot3(INDEI,OMEGA(omega)*ones(1,numel(INDEI)),mean(GEIT{omega}(ensemble+1:2*ensemble,:),1),'b');
    
    subplot(nr,nc,2*nc+5:2*nc+6);hold on;
    plot3(INDIE,OMEGA(omega)*ones(1,numel(INDIE)),mean(GIET{omega}(1:ensemble,:),1),'r');
    plot3(INDIE,OMEGA(omega)*ones(1,numel(INDIE)),mean(GIET{omega}(ensemble+1:2*ensemble,:),1),'b');
end
subplot(nr,nc,2*nc+1:2*nc+2);hold on;
plot3(INDEE(11)*ones(1,numel(OMEGA)),OMEGA,0*ones(1,numel(OMEGA)),'--k')
xlim([INDEE(1) INDEE(end)])

subplot(nr,nc,2*nc+3:2*nc+4);hold on;
plot3(INDEI(11)*ones(1,numel(OMEGA)),OMEGA,0*ones(1,numel(OMEGA)),'--k')
xlim([INDEI(1) INDEI(end)])

subplot(nr,nc,2*nc+5:2*nc+6);hold on;
plot3(INDIE(11)*ones(1,numel(OMEGA)),OMEGA,0*ones(1,numel(OMEGA)),'--k')
xlim([INDIE(1) INDIE(end)])

