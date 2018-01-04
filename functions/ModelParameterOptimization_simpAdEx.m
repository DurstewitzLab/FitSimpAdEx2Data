function [erg,add_para,value,ferr]=ModelParameterOptimization_simpAdEx(I,f,f1,Vnull,Inull,para,para_fix,per)
% performs the actual model optimization

% INPUT:
% * I:              vector containing constant current steps in pA
% * f:              vector containing steady-state firing rates in Hz
% * f1:            	vector containing onset firing rates in Hz
% * Vnull:          vector containing the averaged voltage responses to
%                   step currents below the rhoebase (in mV)
% * Inull:          vector containing the step currents below the rheobase
%                   in pA
% * para:           vector containing the initial values and the lower
%                   boundary for the parameter b
% * para_fix:       vector containing the model parameter numbers that
%                   are kept constant in the final optimization step 
% * per:            tolerated deviation from the estimates of the initial
%                   values

% OUTPUT:
% * erg:            fimal optimzed model parameters 
% * add_para:       parameter values that were kept constant during the
%                   final optimization
% * value:          the initial values (can deviate from the estimated initial 
%                   values) used in the iteration steps (the repititions)
% * ferr:           overall weigthed fitting error in the final
%                   optimization


    % initial values
    Cm=para(1);
    gL=para(2);
    tau=Cm./gL;
    EL=para(3);
    sf=2.0;
    Vr=para(4);
    Vth=para(5);
    bmin=para(6);

    % rename data
    exp_data.I=I;
    exp_data.f=f;

    % variation of initial values
    gL=gL*(1-per)+2*gL*per*rand(1);
    EL=EL*(1-per)+2*EL*per*rand(1);
    sf=sf*(1-per)+2*sf*per*rand(1);
    Vr=Vr*(1-per)+2*Vr*per*rand(1);
    Vth=Vth*(1-per)+2*Vth*per*rand(1);

    % intermediate storage
    value(1)=Cm;
    value(2)=gL;         
    value(3)=EL;   
    value(4)=sf;         
    value(7)=0; 
    value(9)=Vr;
    value(10)=Vth;


    % first step: onset of fI-curve
    OPTIONS=optimset('TolFun',1e-1,'TolX',1e-1,'Display','off');
    [res,ferr]=fminsearch(@null_fit,value([2,3,4,10]),OPTIONS,exp_data.I(1),Vnull([1,end-1,end]),Inull([1,end-1,end]),value([1,5,6,7,8,9]),[1,5,6,7,8,9]);
    value([2,3,4,10])=res;
    value(1)=tau*value(2);
    value(5)=10*value(4)-40;
    if value(9)>value(10)
        value(9)=value(10)-1.0;
    end
    disp('First step: onset and IV');


    % find initial values such that the slope fits roughly
    OPTIONS=optimset('TolFun',1e-1,'TolX',1e-1,'Display','off');
    NumParOp=[6,8,9];
    NumParNOp=[1,2,3,4,5,7,10];
    b_test=[10 50 100 300];
    tcw_test=500;
    for i=1:length(b_test)
        value(8)=b_test(i);
        for j=1:length(tcw_test)
            value(6)=tcw_test(j);
            [erg,ferr]=fminsearch(@slope_fit_both,value(NumParOp),OPTIONS,exp_data.I([1,2,end]),exp_data.f([1,2,end]),f1([1,2,end]),value(NumParNOp),NumParNOp,bmin);
            mtx(i,:)=[erg ferr];
        end
    end
    value(NumParOp)=mtx(min(find(min(mtx(:,length(mtx(1,:))))==mtx(:,length(mtx(1,:))))),1:end-1);
    disp('Second step: slope');


    % final optimization
    OPTIONS=optimset('TolFun',1e-5,'TolX',1e-5,'Display','off');
    u=1; v=1;
    for i=1:10
      if find(i==para_fix)
          add_para(u)=value(i);
          u=u+1;
      else
          para0(v)=value(i);
          v=v+1;
      end;
    end;
    [erg,ferr]=fminsearch(@fit_all3,para0,OPTIONS,exp_data.I,exp_data.f,f1,Vnull,Inull,add_para,para_fix,bmin);
    disp('Final step: optimize the parameters => fit of the onset and steady-sate fI curve as well as the IV curve');

end


% ################### functions ######################

function q=null_fit(value0,I,Vnull,Inull,value,fix)

    u=1; v=1;
    for i=1:10
        if find(i==fix)
            Werte(i)=value(u);
            u=u+1;
        else
            Werte(i)=value0(v);
            v=v+1;
        end;
    end;

    [C,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]=names(Werte);
    I0sim=gL*(Vth-EL)-gL*sf;

    Isim=gL.*(Vnull-EL) - gL*sf*exp((Vnull-Vth)./sf);
    q=sum((Inull-Isim).^2)+(I(1)-I0sim).^2;

end


function q=slope_fit_both(value0,I,fr,fr1,value,fix,bmin)

    u=1; v=1;
    for i=1:10
        if find(i==fix)
            Werte(i)=value(u);
            u=u+1;
        else
            Werte(i)=value0(v);
            v=v+1;
        end;
    end;

    N=length(I);

    % calculate firing-rate
    [C,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]=names(Werte);
    for i=1:N
        re1(i)=FRsimpAdEx(Werte,I(i),b,[],bmin);
        re(i)=FRsimpAdEx(Werte,I(i),[],[],bmin);
    end;

    q=4*sum((re-fr).^2)+sum((re1-fr1).^2);

end


function q=fit_all3(value0,I,fr,fr1,Vnull,Inull,value,fix,bmin)

    u=1; v=1;
    for i=1:10
        if find(i==fix)
            Werte(i)=value(u);
            u=u+1;
        else
            Werte(i)=value0(v);
            v=v+1;
        end;
    end;

    N=length(I);

    % calculate firing-rate
    [C,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]=names(Werte);
    for i=1:N
        re1(i)=FRsimpAdEx(Werte,I(i),b,[],bmin);
        re(i)=FRsimpAdEx(Werte,I(i),[],[],bmin);
    end;

    % calculate IV
    Isim=gL*(Vnull-EL) - gL*sf*exp((Vnull-Vth )./sf);

    % calculate Onset
    I0sim=gL*(Vth-EL)-gL*sf;

    % calculate squared deviation
    q=5*sum((re-fr).^2)+sum((re1-fr1).^2)+4*sum((Inull-Isim).^2)+4*(I(1)-I0sim).^2;

end


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
