function [V_min,V_max,V_th,tau_sp,Vmin,Tkmin,Vmax,Tkmax]=spike_analysis(T,V,T_bounds,plotflag,V_half)
% extracts features of a time series of membrane potentials which contains
% spikes

% INPUT:
% * T:              time vector
% * V:              time series of membrane potential values
% * T_bounds:       vector of time constraints
% * plot_flag:      (1=plot, 0=do not plot)
% * V_half:         value that determines a lower bound for the voltage to
%                   extract the averaged spike peak (e.g. 0)

% OUTPUT:
% * V_min/V_max:     averaged minimum and maximum voltage of a spike
% * V_th:            estimated spiking threshold 
% * tau_sp:          averaged spike width
% * Vmin/Vmax:       time series of all spike peaks and minima
% * Tmin/Tmax:       vector of all time points where spike peaks and minima
%                    are registered

Vmin=[];
Vth=[];
Tkmin=[];
Tth=[];
Vmax=[];
Tkmax=[];
T_up=[];
T_down=[];
tausp=[];

V_act=V(T>T_bounds(1) & T<=T_bounds(2));
T_act=T(T>T_bounds(1) & T<=T_bounds(2));
N=length(V_act);

%Vmax, Tkmax
for k=2:(N-1)
    if (V_act(k)>=V_act(k-1) && V_act(k)>V_act(k+1) && V_act(k) > V_half)
        Vmax=[Vmax,V_act(k)];
        Tkmax=[Tkmax,T_act(k)];
    end;    
end;

% Vmin, Tkmin
for k=2:numel(Tkmax)-2
    q=find(T_act==Tkmax(k));
    p=find(T_act==Tkmax(k+1));
    Vmin(k-1)=min(V_act(q:p));
    r=min(find(V_act(q:p)==Vmin(k-1)));
    Tkmin(k-1)=T_act(r+q-1);
end;

% Vth
dd=0;
for k=1:numel(Tkmax)-3
    q=find(T_act==Tkmin(k));
    while (dd <= 1)
        dd=abs(V_act(q+1)-V_act(q));
        q=q+1;
    end;
    Vth=[Vth,V_act(q-1)];
    Tth=[Tth,T_act(q-1)];
    dd=0;
end;

% average
V_th=mean(Vth(1:floor(length(Vth)/5)));
V_min=mean(Vmin(1:floor(length(Vmin)/5)));
V_max=mean(Vmax(1:floor(length(Vmax)/5)));
if isempty(V_half)
    V_half=(V_min+V_max)/2.0;
end;

for k=1:(N-1)
    if (V_act(k)<=V_half && V_act(k+1)>V_half)
        Tup=((V_half-V_act(k))/(V_act(k+1)-V_act(k)))*(T_act(k+1)-T_act(k))+T_act(k);
        T_up=[T_up,Tup];
    elseif (V_act(k)>V_half && V_act(k+1)<V_half)
        Tdown=((V_half-V_act(k))/(V_act(k+1)-V_act(k)))*(T_act(k+1)-T_act(k))+T_act(k);
        T_down=[T_down,Tdown];
    end;
    if (isempty(T_down)==0 && isempty(T_up)==0)
        if(T_down(length(T_down))-T_up(length(T_up)) > 0)
            diff=T_down(length(T_down))-T_up(length(T_up));
            tausp=[tausp,diff];
        end;
    end;
end;

tau_sp=mean(tausp);

if plotflag
    figure;
    plot(T(T>T_bounds(1) & T<=T_bounds(2)),V_act,'.-b');
    hold on;
    plot(Tkmin,Vmin,'ro');
    plot(Tkmax,Vmax,'go');
    plot(T_up,V_half,'ks');
    plot(T_down,V_half,'ks');
    plot(Tth,Vth,'ys');
end;


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
