clear;
clc;
%% The code is Written by Aref Pariz
% Before using this code, you should have run the data_reader.m. 
% The code will use the data in .MAT file saved in "folder" to find the 
% neurons with smaller and larger MTC with respect to initial values as
% taum0. Then it will find the synaptic weight among those coupled neurons
% and will find the distribution of synaptic weights. 



N=10000;
Ne=0.8*N;
Ni=N-Ne;
OMEGA=[10 15 20 25]; % Stimulation frequencies.
AMP=10; % Stimulation amplitude(s).
folder='data/'; % Data are kept in this folder
ca=1; % This value defines the name of the file that data will be save into.

PART=21:20; % Post-stimulation epoch. This value is being used to take average over time from synaptic weights
EE=1;EI=1;IE=1; % Synaptic weights that have been modified. same as STDP values in Fortran code.

gee=0.00005; % Initial synaptic weights between E -> E neurons
gei=0.00005; % Initial synaptic weights between E -> I neurons
gie=0.00025; % Initial synaptic weights between I -> E neurons
gii=0.00025; % Initial synaptic weights between I -> I neurons



INDEE=linspace(0,2*gee,21); 
INDEI=linspace(0,2*gei,21);
INDIE=linspace(0,2*gie,21);

ensemble=5;

% In case in the program the stimulation amplitude is changing, here the
% FOR loop is over amp and the OMEGA is the stimulation amplitudes. If the
% amplitude of stimulation is changing, the "amp" should be changed to
% "omega" and the AMP should contain the stimulation amplitudes and ca=2.

for amp=AMP % Stimulation amplitude
    nomega=0;
    
    if EE
        GEET=cell(1,numel(OMEGA));
    end
    if EI
        GEIT=cell(1,numel(OMEGA));
    end
    if IE
        GIET=cell(1,numel(OMEGA));
    end
    for omega=OMEGA
        disp(omega);
        nomega=nomega+1;
        if EE
            valxEE_less=zeros(ensemble,numel(INDEE));
            valxEE_more=zeros(ensemble,numel(INDEE));
        end
        if EI
            valxEI_less=zeros(ensemble,numel(INDEI));
            valxEI_more=zeros(ensemble,numel(INDEI));
        end
        if IE
            valxIE_less=zeros(ensemble,numel(INDIE));
            valxIE_more=zeros(ensemble,numel(INDIE));
        end
        for ens=1:ensemble
            fname0=[folder,'params_omega',num2str(omega),'_amp',num2str(amp),'_ens',num2str(ens),'.mat'];
            if exist(fname0,'file')
                load(fname0,'AB','taum');
                
                taum0=10;
                dtaum=2;
                taumE=taum(1:Ne);
                taumI=taum(Ne+1:N);
                
                taumE_less=find(taumE<taum0-dtaum);
                taumE_more=find(taumE>taum0+dtaum);
                
                taumI_less=Ne+find(taumI<taum0-dtaum);
                taumI_more=Ne+find(taumI>taum0+dtaum);
                
                taumE_initial = find(abs(taumE-taum0)<0.5);
                taumI_initial = Ne+find(abs(taumI-taum0)<0.5);
                
                syn_numE=numel(taumE_initial);
                syn_numI=numel(taumI_initial);
                %%
                JsumEE_less=zeros(syn_numE,numel(taumE_less));
                JsumEE_more=zeros(syn_numE,numel(taumE_more));
                
                JsumEI_less=zeros(syn_numE,numel(taumI_less));
                JsumEI_more=zeros(syn_numE,numel(taumI_more));
                
                JsumIE_less=zeros(syn_numI,numel(taumE_less));
                JsumIE_more=zeros(syn_numI,numel(taumE_more));
                for part=PART
                    fname1=[folder,'result_omega',num2str(omega),...
                        '_amp',num2str(amp),'_ens',num2str(ens),...
                        '_part',num2str(part),'.mat'];
                    load(fname1,'Jsave');
                    J_regenerate
                    
                    if EE
                        n1=0;
                        for ii=taumE_initial'
                            n1=n1+1;
                            for jj=1:numel(taumE_less)
                                pre=ii;
                                post=taumE_less(jj);
                                if A(post,pre)
                                    JsumEE_less(n1,jj) = JsumEE_less(n1,jj) + J(post,pre);
                                end
                            end
                        end
                        
                        n1=0;
                        for ii=taumE_initial'
                            n1=n1+1;
                            for jj=1:numel(taumE_more)
                                pre=ii;
                                post=taumE_more(jj);
                                if A(post,pre)
                                    JsumEE_more(n1,jj) = JsumEE_more(n1,jj) + J(post,pre);
                                end
                            end
                        end
                    end
                    if EI
                        n1=0;
                        for ii=taumE_initial'
                            n1=n1+1;
                            for jj=1:numel(taumI_less)
                                pre=ii;
                                post=taumI_less(jj);
                                if A(post,pre)
                                    JsumEI_less(n1,jj) = JsumEI_less(n1,jj) + J(post,pre);
                                end
                            end
                        end
                        
                        n1=0;
                        for ii=taumE_initial'
                            n1=n1+1;
                            for jj=1:numel(taumI_more)
                                pre=ii;
                                post=taumI_more(jj);
                                if A(post,pre)
                                    JsumEI_more(n1,jj) = JsumEI_more(n1,jj) + J(post,pre);
                                end
                            end
                        end
                        
                    end
                    if IE
                        n1=0;
                        for ii=taumI_initial'
                            n1=n1+1;
                            for jj=1:numel(taumE_less)
                                pre=ii;
                                post=taumE_less(jj);
                                if A(post,pre)
                                    JsumIE_less(n1,jj) = JsumIE_less(n1,jj) + J(post,pre);
                                end
                            end
                        end
                        
                        n1=0;
                        for ii=taumI_initial'
                            n1=n1+1;
                            for jj=1:numel(taumE_more)
                                pre=ii;
                                post=taumE_more(jj);
                                if A(post,pre)
                                    JsumIE_more(n1,jj) = JsumIE_more(n1,jj) + J(post,pre);
                                end
                            end
                        end
                        
                    end
                end
                if EE
                    JsumEE_less=JsumEE_less./numel(PART);
                    JsumEE_more=JsumEE_more./numel(PART);
                    
                    xEE=JsumEE_less(:);
                    xEE(xEE==0)=[];
                    [aux,~]=hist(xEE,INDEE);
                    valxEE_less(ens,:)=aux;
                    
                    xEE=JsumEE_more(:);
                    xEE(xEE==0)=[];
                    [aux,~]=hist(xEE,INDEE);
                    valxEE_more(ens,:)=aux;
                end
                if EI
                    JsumEI_less=JsumEI_less./numel(PART);
                    JsumEI_more=JsumEI_more./numel(PART);
                    
                    xEI=JsumEI_less(:);
                    xEI(xEI==0)=[];
                    [aux,~]=hist(xEI,INDEI);
                    valxEI_less(ens,:)=aux;
                    
                    xEI=JsumEI_more(:);
                    xEI(xEI==0)=[];
                    [aux,~]=hist(xEI,INDEI);
                    valxEI_more(ens,:)=aux;
                end
                if IE
                    JsumIE_less=JsumIE_less./numel(PART);
                    JsumIE_more=JsumIE_more./numel(PART);
                    xIE=JsumIE_less(:);
                    xIE(xIE==0)=[];
                    [aux,~]=hist(xIE,INDIE);
                    valxIE_less(ens,:)=aux;
                    
                    xIE=JsumIE_more(:);
                    xIE(xIE==0)=[];
                    [aux,~]=hist(xIE,INDIE);
                    valxIE_more(ens,:)=aux;
                end
            else
                warning(['file for omega=',num2str(omega),' does not exist'])
            end
        end
        if EE
            GEET{1,nomega}=[valxEE_less;valxEE_more];
        end
        if EI
            GEIT{1,nomega}=[valxEI_less;valxEI_more];
        end
        if IE
            GIET{1,nomega}=[valxIE_less;valxIE_more];
        end
    end
    %%
    clear A Jsave J AB
    switch ca
        case 1
            fname=['result_amp',num2str(amp),'.mat'];
            save(fname,'GEIT','GEET','GIET','amp','OMEGA','PART','INDEE','INDEI','INDIE',...
                'gee','gei','gie','folder','ensemble')
        case 2
            fname=['result_omega',num2str(omega),'.mat'];
            save(fname,'GEIT','GEET','GIET','amp','OMEGA','PART','INDEE','INDEI','INDIE',...
                'gee','gei','gie','folder','ensemble')
    end
end
%%
