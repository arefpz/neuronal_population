%% Written by Aref Pariz
% simulate the activity of coupled neurons in presence of stimulation with
% frequency omega and neurons' MTC as "taum"
clc
clear
rng('shuffle');
% FOr older version, uncomment below
% stream = RandStream('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(stream);

%% Constants
sec=100; % in second
dt=0.1; % in ms
range=1000*sec/dt; % TIME STEPS
N=2; % number of neurons per network

% connetions and weights
Nep=1; % number of excitatory input for each neuron
Nip=1-Nep; % number of inhibitory input for each neuron
Ne=Nep*N;   Ni=N-Ne;
taum(1,1)=10;
taum(2,1)=20;
v_thr=-54*ones(N,1); % Threshold for spike
I0=5.5*ones(N,1);
A=[0,1;1,0];
tref=2/dt;
rho=zeros(N,range,'logical');
% Iext=external_current(numnet*N,range,dh,r,taudIext);
sqrt_dt=1/sqrt(dt);
sigma=1;

V=zeros(N,range);
vrest=-60*ones(N,1);

Ap=0.02;Am=0.01;
tauP=10;tauM=10;
DC = 0;
S=zeros(N,1);
taud=3;
gmin=0.01;gmax=0.2;
ensemble=10;
%% main
folder='data/';
try
    cd data/
    cd ..
catch
    mkdir data
end

n1=0;
OMEGA=0:2:50;
TAUM=4:1:20;
Gend12=zeros(numel(TAUM),numel(OMEGA));
Gend21=Gend12;
amp=0.5; % 1
for omega=OMEGA
    n1=n1+1;n2=0;
    display(omega)
    for tau=TAUM
        taum(2,1)=tau;
        n2=n2+1;
        display(amp)
        Isignal=amp * sin(2*pi*omega*(0:range-1)/10000);
        for ens=1:ensemble
            fname=[folder,'result_omega',num2str(omega),'_taum',num2str(tau),'_ens',num2str(ens),'.mat'];
            if ~exist(fname,'file')
                v=zeros(N,1);
                tspike=zeros(N,1);
                G12=zeros(1,range);
                G21=zeros(1,range);
                g21=0.1;g12=0.1;
                g0=g12;
                RI=zeros(N,1);
                Inoise=sqrt_dt*normrnd(0,sigma,N,range);
                t_ref=zeros(2,1);
                for ii=1:range
                    v = v + dt * (vrest-v + I0 + DC +Inoise(:,ii) + Isignal(ii)+RI) ./taum;
                    fired=(v>=v_thr); % if v is above threshod and previous v is less than threshod, neuron fired
                    t_ref(fired)=ii+tref;
                    O1=t_ref>=ii;
                    v(O1)=vrest(O1);
                    rho(:,ii)=fired;
                    tspike(fired)=ii;
                    S(fired) =S(fired) +1;
                    S=S+dt*(-S/taud);
                    if fired(1)
                        dT=dt*(tspike(1)-tspike(2));
                        g21=g21+Ap*(1-g21/gmax)*exp(-dT/tauP);
                        g12=g12-Am*(g12/g0)*exp(-dT/tauM);
                    elseif fired(2)
                        dT=dt*(tspike(2)-tspike(1));
                        g21=g21-Am*(g21/g0)*exp(-dT/tauM);
                        g12=g12+Ap*(1-g12/gmax)*exp(-dT/tauP);
                    end
                    if g12<gmin
                        g12=gmin;
                    elseif g12>gmax
                        g12=gmax;
                    end
                    if g21<gmin
                        g21=gmin;
                    elseif g21>gmax
                        g21=gmax;
                    end
                    G21(ii)=g21;
                    G12(ii)=g12;
                    V(:,ii)=v;
                    RI = [S(2) * g21;S(1) * g12];
%                     Isyn(:,ii)=RI;
%                     SS(:,ii)=S;
                end
                save(fname,'-v7.3','G12','G21','V','rho')
            end
        end
    end
end

%%
