function exp_data=define_onset(I,f,f1)
% extracts the onset of the f-I curves from the real data and provides data extended by the estimated onset 

% INPUT:
% * I:          vector containing the applied currents in pA
% * f:          vector containing the steady-state firing rates in Hz
% * f1:         vector containing the onset firing rates in Hz

% OUTPUT:
% * exp_data:   structure with "new" I, f and f1 (extended by the onest)


    % define relevant experimental data:
    l=min(find(f~=0));
    N=length(I);
    j=0;
    if(l==1)
        s=l;
    else
        s=l-1;
    end;
    for i=s:N
        j=j+1;
        exp_data.I(j)=I(i);
        exp_data.f(j)=f(i);
        exp_data.f1(j)=f1(i); 
    end;


    % fit log function
    if(l~=1)
        var(2)=-exp_data.I(1);
    else
        res=polyfit(exp_data.I(1:2),exp_data.f(1:2),1);
        var(2)=res(2)./res(1);
    end;
    var(1)=(exp_data.f(3)-exp_data.f(2))/(log(exp_data.I(3)+var(2))-log(exp_data.I(2)+var(2)));
    var(3)=exp_data.f(2)-var(1)*log(exp_data.I(2)+var(2));
    OPTIONS=optimset('TolFun',1e-5,'TolX',1e-5);
    if(l==1)
        [res,fval]=fminsearch(@fit_log,var,OPTIONS,exp_data.I(1:end),exp_data.f(1:end));
    else
        [res,fval]=fminsearch(@fit_log,var,OPTIONS,exp_data.I(2:end),exp_data.f(2:end));
    end;
    I0=exp(-res(3)./res(1))-res(2);


    % correction if estimated onset is greater than the I-value
    % corresponding to the first non-zero firing rate
    if(I0 > exp_data.I(2))
        if(l~=1)
            I0=exp_data.I(1);
        else
            %linear fit
            m=(exp_data.f(1)-exp_data.f(2))./(exp_data.I(1)-exp_data.I(2));
            n=exp_data.f(1)-m*exp_data.I(1);
            I0=-n/m;
        end;
    end;

    % correction if estimated onset is smaller than the I-value corresponding
    % to the last zero firing rate
    if(l~=1)
        if I0<exp_data.I(1)
            I0=exp_data.I(1);
        end
    end

    % display estimated onset
    str=['estimated onset: ' num2str(I0) ' pA'];
    disp(str);

    % if necessary, extend the I (and f) vector(s) by the onset
    if(l==1)
        exp_data.I=[I0,exp_data.I(1:end)];
        exp_data.f=[0,exp_data.f(1:end)];
        exp_data.f1=[0,exp_data.f1(1:end)];
    else
        exp_data.I=[I0,exp_data.I(2:end)];
    end

end


% ################### function/s ######################


function q=fit_log(a,xd,yd)
    % simple log model for the f-I curve

    y=a(1)*log(xd+a(2))+a(3);
    q=sum((y-yd).^2);
end


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
