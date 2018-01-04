function SimPar = set_parameters
% sets all parameters needed (e.g. to run the model simulation
% voltage_AnalyticalModel_simpAdEx.c)

% ------------------------ Set meta data ----------------------------------
SimPar.name = 'all_parameters';
SimPar.comment = 'none';

% ------------------------- Set neuron parameters -------------------------
CellPar = zeros(11,1);         % (mean) neuron parameters

CellPar(1,:)  = 0;             % I: external driving current (NOT used!)

CellPar(2,:)  = 200;           % Cm: membrane capacity (setting membrane time constant)
CellPar(3,:)  =  10;           % gL: leak current peak conductance
CellPar(4,:)  = -70;           % EL: leak current reversal potential
CellPar(5,:)  =   2.0;         % sf: slope factor
CellPar(6,:)  = -10.0;         % Vup: upswing reset threshold
CellPar(7,:)  =  30.0;         % tcw: tau_w
CellPar(8,:)  =   0.0;         % a: bifurcation parameter
CellPar(9,:)  =   2.0;         % b: bifurcation parameter
CellPar(10,:) = -58;           % Vr: reset potential
CellPar(11,:) = -50;           % Vthr: Threshold for spiking

% ----------------------- Pass parameters to SimPar -----------------------
SimPar.CellPar = CellPar;

% -------------------------- Save SimPar -------------------------------
SimFile=[SimPar.name '.mat'];
save(SimFile,'SimPar');


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
