function [res,index1,index2,b] = GenerateTrainingSet(folder,FileIn,FileOut,T_bounds_dep,T_bounds_hyper,thrs,fI_idx_all)
% extracts the steady-state and onset f-I curves as well as the
% sub-rheobase I-V curve from experimental data

% INPUT:
% * folder:             directory containing the raw data set 
%                       (responses to constant step currents)
% * FileIn:             file name of the raw data set
% * FileOut:            file name of the training set 
% * T_bounds_dep:       time constraints for averaging 
%                       (responses to depolarizing current steps)
% * T_bounds_hyper:     time constraints for averaging 
%                       (responses to hyperpolarizing current steps)
% * thres:              spiking threshold in mV
% * fI_idx_all:         vector containing all indices of fI curves 

% OUTPUT:
% * res:                structure with current steps, steady-state and onset firing
%                       rates as well as the IV relation (the training set) 
% * index1:             indices (row numbers) containing the
%                       hyperpolarizing steps (one index = 4 hyperpol. steps)
% * index2:             indices (row numbers) containing the depolarizing
%                       steps (to compute the f-I curves)
% * b:                  lower boundary for the model parameter b

% IMPORTANT: It is assumed that I is provided in pA (steps from 25 pA up to
% 700 pA in steps of 25 pA) and 4 hyperpolarizing steps were recorded 
% (-200, -150, -100 and -50 pA). Furthermore, the time is presumed to be 
% recorded in sec and the voltage resonse in V! Please make sure that your 
% data meet the criteria or adjust the implementation below.

    warning off;

    % load file
    load([folder '/d_' FileIn]);

    % find relevant indices and rename them
    index1=fI_idx_all-1;    % indices of I-V curve data - assumed to be right before each f-I curve
    index2=fI_idx_all;      % indices of f-I curve data

    % steady-state firing rates
    thr=thrs*ones(1,length(index2));
    I=25:25:700;
    for j=1:length(index2) % loop over all f-I curves
        fI_idx = index2(j);
        res(j).I = I(1:length(data.sweep(fI_idx,data.sweep(fI_idx,:)>0)));
        res(j).f = zeros(length(data.sweep(fI_idx,data.sweep(fI_idx,:)>0)),2);
        for i=1:length(data.sweep(fI_idx,data.sweep(fI_idx,:)>0))
            res_act = fI_calc2(data.t{fI_idx,i}, res(j).I(i), data.v{fI_idx,i}*1000, T_bounds_dep, thr(j));
            res(j).f(i,1:2) = res_act.f;
        end;
    end

    % onset firing rates and IV relation
    for i=1:length(index2) % loop over all f-I curves
        res(1,i).Inull=[];
        Nm=min(find(data.sweep(index2(i),:)==0))-1; % number of data points of respective f-I curve
        if isempty(Nm)
            Nm=length(data.sweep(index2(i),:));
        end
        N0=min(find(res(1,i).f(:,1)~=0)); % first non-zero firing rate
        B=0;
        for j=1:Nm
            if j<N0
                res(1,i).Inull(j)=res(1,i).I(j);
                res(1,i).Vnull(j)=mean(data.v{index2(i),j}(data.t{index2(i),j}>T_bounds_dep(1) & data.t{index2(i),j}<T_bounds_dep(2))*1000);
                res(1,i).f1(j)=0;
            else
                train=sp_train(data.t{index2(i),j},data.v{index2(i),j}*1000,thrs,0);
                if length(train)>=5
                    ISI=[train(1) diff(train)'];
                    f=1./ISI';
                    r=polyfit(train(2:3),f(2:3),1);
                    x0=r(2)-res(1,i).f(j,1);
                    p=fminsearch(@func_exp,[x0 r(1)./x0 res(1,i).f(j,1)],[],train(3:end),f(3:end));
                    res(1,i).f1(j)=p(1)*exp(p(2)*train(1))+p(3);
                    B=B+abs(p(2));
                else
                    res(1,i).f1(j)=0;
                end
            end
        end
        b(i)=ceil(B./(length(res(1,i).f1)));
    end

    % IV relation: 4 hyperpolarizing steps
    for i=1:length(res) % loop over all I-V curves
        for j=1:4
            V0(j)=mean(data.v{index1(i),j}(data.t{index1(i),j}>T_bounds_hyper(1) & data.t{index1(i),j}<T_bounds_hyper(2))*1000);
        end
        if isempty(res(1,i).Inull)
            res(1,i).Inull=[-200 -150 -100 -50];
            res(1,i).Vnull=[V0];
        else
            res(1,i).Inull=[-200 -150 -100 -50 res(1,i).Inull];
            res(1,i).Vnull=[V0 res(1,i).Vnull];
        end
    end

    % save
    save([folder '/' FileOut '.mat'],'res','index1','index2','b');
    disp('Training set is computed and saved');
end

% ################### functions ######################

function q = func_exp(x0,xd,yd)

    y=x0(1).*exp(x0(2).*xd)+x0(3);
    q=sum((y-yd).^2);

end


function res = fI_calc2(T, I, V, T_bounds, Vth)
% calculates f-I curves from data consisting of voltages traces V and
% injected currents I. Format of both matrices must be M(time,trials).
% I_on and I_off denote the time window within which the current injection
% is considered (2.5 and 27.5 sec, by default). T is the time vector.

    N = length(I(1,:));         % number of trials with different I
    I_on = 2.5;
    I_off = 27.5;

    for i=1:N
        if length(I(:,i))>1
            res.I(i,1) = mean(I(T>=I_on & T<I_off,i));
            res.I(i,2) = std(I(T>=I_on & T<I_off,i));
        else
            res.I(i,1) = I;
            res.I(i,2) = 0;
        end

        k=find(T>T_bounds(1) & T<=T_bounds(2));
        ISI_act=IntSpInv2(T(k),V(k,i),Vth,0);

        if isempty(ISI_act)
            res.f(i,1:2) = 0;
        else
            res.f(i,1) = mean(1./ISI_act);
            res.f(i,2) = std(1./ISI_act);
        end
    end

end


function ISI = IntSpInv2(T,V,Vth,plot_flag)

    k=find(V(2:end)>=Vth & V(1:end-1)<Vth);
    Tsp=T(k)+(T(k+1)-T(k)).*(Vth-V(k))./(V(k+1)-V(k));
    ISI=diff(Tsp);
    if plot_flag
        plot(T,V,'b');
        hold on, plot(Tsp,Vth,'ro');
    end

end


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg
% and BCCN Heidelberg-Mannheim