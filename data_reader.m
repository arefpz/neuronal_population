%% Reading data from .txt file generate by Fortran
clc
clear
filename = 'variable.txt';
formatSpec = '%13f%17f%6s%2s%2s%2s%2s%s%[^\n\r]';
folder='';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);

dataArray{3} = strtrim(dataArray{3});
dataArray{4} = strtrim(dataArray{4});
dataArray{5} = strtrim(dataArray{5});
dataArray{6} = strtrim(dataArray{6});
dataArray{7} = strtrim(dataArray{7});
dataArray{8} = strtrim(dataArray{8});
fclose(fileID);
dataArray([1, 2]) = cellfun(@(x) num2cell(x), dataArray([1, 2]), 'UniformOutput', false);
var = [dataArray{1:end-1}];

%%
PART=30; % number of times data has been written
ensemble=5; % number of ensemble
N=10000; % Number of neurons
each_part=10000; % Each part time steps
for ens=1:ensemble
    part=0;
    disp(ens)
    log=readtable('log.txt');
    selected_neurons=single(load(['selcN_part',num2str(part),'_ens',num2str(ens),'.txt']));
    selected_synapses=single(load(['selcS_part',num2str(part),'_ens',num2str(ens),'.txt']));
    taum=single(load(['taum_part',num2str(part),'_ens',num2str(ens),'.txt']));
    J0=single(load(['ABsynapse0_part',num2str(part),'_ens',num2str(ens),'.txt']));
    AB=single(load(['ABval_part',num2str(part),'_ens',num2str(ens),'.txt']));
    vthr=single(load(['vthreshold_part',num2str(part),'_ens',num2str(ens),'.txt']));
    fname0=['params_omega',num2str(var{1}),'_amp',num2str(var{2}),'_ens',num2str(ens),'.mat'];
    save(fname0,'-v7.3','J0','log','vthr','selected_neurons','selected_synapses', 'taum', 'AB');
    for part=1:PART
        disp(part)
        rho=zeros(N,each_part,'logical');
        fname1=['rho_part',num2str(part),'_ens',num2str(ens),'.txt'];
        rhox=load(fname1);
        for jj=1:N
            for ii=1:size(rhox,2)
                if rhox(jj,ii) ~= 0
                    rho(jj,rhox(jj,ii)) = 1;
                end
            end
        end
        
        
        Ve=single(load([folder,'VEneu_part',num2str(part),'_ens',num2str(ens),'.txt']));
        Vi=single(load([folder,'VIneu_part',num2str(part),'_ens',num2str(ens),'.txt']));
        
        fnameGee=[folder,'Gsumee_part',num2str(part),'_ens',num2str(ens),'.txt'];
        Gsumee=single(load(fnameGee));
        
        fnameGei=[folder,'Gsumei_part',num2str(part),'_ens',num2str(ens),'.txt'];
        Gsumei=single(load(fnameGei));
        
        fnameGie=[folder,'Gsumie_part',num2str(part),'_ens',num2str(ens),'.txt'];
        Gsumie=single(load(fnameGie));
        
        
        
        Jsave=single(single(load([folder,'ABsynapse_part',num2str(part),'_ens',num2str(ens),'.txt'])));
        
        V=single(load([folder,'Vsave_part',num2str(part),'_ens',num2str(ens),'.txt']));
        stdp=single(load([folder,'stdpc_part',num2str(part),'_ens',num2str(ens),'.txt']));
        rho = logical(rho);
        
        fname1=['result_omega',num2str(var{1}),'_amp',num2str(var{2}),'_ens',num2str(ens),'_part',num2str(part),'.mat'];
        
        save(fname1,'-v7.3','Ve','Vi','Jsave','V','Gsumee','Gsumei','Gsumie','stdp','rho');
    end
end
