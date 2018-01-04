function train=sp_train(T,V,Vth,plot_flag)
% extracts a spike train from a time series of membrane potentials

% INPUT:
% * T:              time vector
% * V:              time series of membrane potential values
% * Vth:            spiking threshold
% * plot_flag:      (1=plot, 0=do not plot)

% OUTPUT:
% * train:          spike train (vector of spike times)

train=[];
k=find(V(2:end)>=Vth & V(1:end-1)<Vth);
Tsp=T(k)+(T(k+1)-T(k)).*(Vth-V(k))./(V(k+1)-V(k));
train=[train,Tsp];

if plot_flag
   figure;
   plot(T,V,'b');
   if ~isempty(Tsp)
     hold on, plot(Tsp,Vth,'ro');
   end;
end;


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
