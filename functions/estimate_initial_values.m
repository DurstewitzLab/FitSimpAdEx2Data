function initial=estimate_initial_values(file1,file2,idx1,idx2,obs,I0,flag)
% estimates initial values from the raw data (for the optimization procedure)

% INPUT:
% * file1:          directory and file name of raw data set 
%                   (responses to constant step currents)
% * file2:          directory and file name of the training set
% * idx1:         	indices for hyperpolarizing steps (I-V curves)
% * idx2:           indices for depolarizing steps (f-I curves)
% * obs:            specifies the number of the trial to be used
% * I0:             onset current in pA
% * flag:           determines whether voltage response is plotted (1=plot, 0=do not plot) 

% OUTPUT:
% * initial:        vector containing the estimates for [C,gL,EL,Vr,Vth,tau]

% IMPORTANT: It is assumed that the depolarizing current steps were applied for at
% least 30 seconds (with the first and last 2.5 sec I = 0 pA) and the 
% hyperpolarizing current steps for 1.1 seconds (with the first and last 0.05 sec I = 0 pA)!
% Furthermore, the time is presumed to be recorded in sec and the voltage 
% resonse in V! Please make sure that your data meet the criteria or adjust the implementation below.

% load files
load(file1)
load(file2)

% part to extract EL,C,gL,tau,Vth
N = 4;
L_index=length(idx1);
fI_idx_all=[];
for i=1:L_index
    fI_idx_all=[fI_idx_all idx1(i)];
end
R_in = zeros(1,length(fI_idx_all));
EL = zeros(1,length(fI_idx_all));
tau = zeros(1,length(fI_idx_all));
tau_res = zeros(N,length(fI_idx_all));
for j=1:length(fI_idx_all)
    R_row = fI_idx_all(j);
    R_res = zeros(N,4);
    R_res(:,1) = [-200 -150 -100 -50]';
    for i=1:N
        % R_in, EL
        R_res(i,2) = (mean(data.v{R_row,i}(data.t{R_row,i}>0.2 & data.t{R_row,i}<=1.0)) - mean(data.v{R_row,i}(data.t{R_row,i}<0.05)))*1000;
        R_res(i,3) = mean(data.v{R_row,i}(data.t{R_row,i}>0.2 & data.t{R_row,i}<=1.0))*1000;
        R_res(i,4) = mean(data.v{R_row,i}(data.t{R_row,i}<0.05))*1000;
        
        % tau
        t = data.t{R_row,i}(data.t{R_row,i}>0.055 & data.t{R_row,i}<=0.1)*1000;
        y = data.v{R_row,i}(data.t{R_row,i}>0.055 & data.t{R_row,i}<=0.1)*1000;
        x0 = [abs(R_res(i,2)) 30 R_res(i,3)]';
        p = fminsearch(@(x) sum(abs(x(1)*exp(-(t-55)/x(2))+x(3)-y)), x0);
        tau_res(i,j) = p(2);
    end;
    
    R_res(:,5) = R_res(:,3) - R_res(:,4) + mean(R_res(:,4));
    
    par = polyfit(R_res(:,1), R_res(:,5), 1);
    R_in(j) = par(1)*1000;      
    EL_all(j) = par(2);   
    Vth_all(j)=par(1)*I0+par(2);
end;
tau_all = mean(tau_res);
C_all = 1000*tau_all./R_in;
gL_all= 1000./R_in;

tau=tau_all(obs);
EL=EL_all(obs);
gL=gL_all(obs);
C=C_all(obs);
Vth=Vth_all(obs);

% part to extract Vr
V_slice=[];
T_slice=[];

fI_row = idx2(obs);
nummin=min(find(res(1,obs).f(:,1)~=0));
nummax=max(find(res(1,obs).f(:,1)~=0));

for i=nummin:nummax
    [V_min,V_max,V_th,tau_sp,Vmin,Tkmin,Vmax,Tkmax]=spike_analysis(data.t{fI_row,i},data.v{fI_row,i}*1000,[2 20],flag,0);

    Vr(i)=V_min;
end;

Vr=mean(Vr(nummin:nummax));
if Vr>Vth
    Vr=Vth-0.5;
end;


% check values and correct if necessary
if (tau<=0 || isnan(tau) || isinf(tau))
    tau=30;
    gL=5;
    C=tau*gL;
end
if (C<=0 || isnan(C) || isinf(C))
    C=200;
    gL=C/tau;
end
if (gL<=0 || isnan(gL) || isinf(gL))
    gL=5;
    C=tau*gL;
end
if (EL>=0 || isnan(EL) || isinf(EL))
    EL=-70;
end
if (Vth<EL || isnan(Vth) || isinf(Vth))
    Vth=-50;
end
if (Vr>=Vth || isnan(Vr) || isinf(Vr))
    Vr=Vth-10;
end


% transfer into output varaible
initial=[C,gL,EL,Vr,Vth,tau];


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
