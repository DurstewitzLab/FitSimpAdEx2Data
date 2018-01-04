function [initial0,start,results]=FitTrainingSet(file,filename,dpt,obs,para_fix,numrep,dev,plot_flag,saveplot_flag,savedata_flag,AddName)
% main function for fitting the parameters of the simpAdEx

% INPUT:
% * file:               directory where the files will be saved
% * filename:           name of the data set
% * dpt:                number of used data points in the f-I curve (needs to be restricted if firing rates decrease at higher currents)
% * obs:                index of the f-I curve to be used for fitting 
% * para_fix:           vector containing the model parameter numbers that are kept constant in the final optimization step 
% * numrep:             number of repetitions of the optimization procedure
% * dev:                tolerated deviation from the estimates of the initial values (e.g. 0.2)
% * plot_flag:          (1=plot fit, 0=do not show fit)
% * saveplot_flag:      (1=save plot ["plot_flag" has to be set to 1 as well], 0=do not save plot) 
% * savedata_flag:      (1=save data, 0=do not save data) 

% OUTPUT:
% * initial0:           initial values of the model parameters, estimated
%                       directly from the data
% * start:              initial values (which may deviate from initial0) used in 
%                       the iteration steps (the number of iteration steps 
%                       (repetitions) is specified in "numrep", see below)
% * results:            matrix containing the optimized model parameters 
%                       (row = iteration number, column = parameters and 
%                       weighted fitting error, in detail: 1=C, 2=gL, 3=EL,
%                       4=sf, 5=Vup, 6=tcw, 7=a (a=0), 8=b, 9=Vr, 10=Vth, 
%                       11=weighted fitting error)


%********************
%**** Load data *****
%********************

dir1=[file 'd_' filename];
load(dir1); 
dir2=[file 'Training_' filename];
load(dir2); 
bmin=b;


%*********************************
%**** Define f-I curve onset *****
%*********************************

exp_data=define_onset(res(1,obs).I(1:dpt),res(1,obs).f(1:dpt,1),res(1,obs).f1(1:dpt));


%**********************************
%**** Estimate initial values *****
%**********************************

initial0=estimate_initial_values(dir1,dir2,index1,index2,obs,exp_data.I(1),0);
disp('estimate initial values');


%***********************
%**** Main program *****
%***********************

if plot_flag
    figure;
    hold on;
end

for i=1:numrep
    
    % show progress
    disp(['%%% Iteration ' num2str(i) ' %%%']);
    
    % paramter optimization
    [erg,add_para,value,ferr]=ModelParameterOptimization_simpAdEx(exp_data.I,exp_data.f,exp_data.f1,res(1,obs).Vnull,res(1,obs).Inull,[initial0(1:end-1) bmin(obs)],para_fix,dev);
    
    % compute the fI-curve with optimized parameters
    SimPar = parameter;
    u=1; v=1;
    for l=1:10
        if find(l==para_fix)
            SimPar.CellPar(l+1,:)=add_para(u);
            u=u+1;
        else
            SimPar.CellPar(l+1,:)=erg(v);
            v=v+1;
        end;
    end;
    
    % transfer results into a matrix
    for j=1:10
        results(i,j)=SimPar.CellPar(j+1,:);
        start(i,j)=value(j);
    end;
    results(i,11)=ferr;

    % compute f-I curve(s) and I-V curve for comparsion --> plot
    if plot_flag==1        
        N=length(res(1,obs).I(1:dpt));
        [C,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]=names(SimPar.CellPar(2:11,:));

        for k=1:N
            re1(k)=FRsimpAdEx(SimPar.CellPar(2:11,:),res(1,obs).I(k),b,[],bmin(obs));
            re(k)=FRsimpAdEx(SimPar.CellPar(2:11,:),res(1,obs).I(k),[],[],bmin(obs));
        end;

        Isim=gL.*(res(1,obs).Vnull-EL) - gL*sf*exp((res(1,obs).Vnull-Vth)./sf);
    
        subplot(numrep,2,2*i-1), plot(res(1,obs).I(1:dpt),re,'ro-',res(1,obs).I(1:dpt),re1,'bo-');
        hold on, plot(res(1,obs).I(1:dpt),res(1,obs).f(1:dpt,1), '.-k',res(1,obs).I(1:dpt),res(1,obs).f1(1:dpt),'.-k');
        title([num2str(i),'. try of intial values']);
        subplot(numrep,2,2*i), plot(res(1,obs).Vnull,res(1,obs).Inull,'.-k',res(1,obs).Vnull,Isim,'.-r');
    end
end
disp('Optimization completed');


%********************************
%**** save data and figures *****
%********************************

% save plot
if saveplot_flag && plot_flag
    if exist('AddName', 'var')==1
        saveas(gcf,[file 'plot_TrainingFit_' filename '_' AddName '_obs' num2str(obs) '_dp' num2str(dpt) '.fig']);
    else
        saveas(gcf,[file 'plot_TrainingFit_' filename '_obs' num2str(obs) '_dp' num2str(dpt) '.fig']);
    end
end

% save data
if savedata_flag
    if exist('AddName', 'var')==1
        save([file 'para_' filename '_' AddName '_obs' num2str(obs) '_dp' num2str(dpt) '.mat'],'start','results','initial0');
    else
        save([file 'para_' filename '_obs' num2str(obs) '_dp' num2str(dpt) '.mat'],'start','results','initial0');
    end
end


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
