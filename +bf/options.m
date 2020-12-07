function machinePref = options(opt)
% Get/Set options associated with the bayesFactor toolbox
% This function is called by bf internal functions. 
% User can fine tune options by first calling this function without 
% input arguments, then modifying the struct that is returned and then
% calling this function again, with the modified struct as its input.
% EG. To use 20 parpool workers by default, use:
% opt= bf.options
% opt.nrWorkers = 20;
% bf.options(opt);
%
% Please note that changes to the options are permanent (i.e. they will be
% used again in the next Matlab session on the same machine).



%% Current factory default parameters for the Monte Carlo integration (see internal.mcIntegrate) and parallel execution
factoryDefault.minG = 0.0001;  % Smalles g value (don't include zero).
factoryDefault.maxG = 10000;   % Largest g value to intergare
factoryDefault.stepG = 0.05;   % Step size for the g integration.
factoryDefault.nrSamples = 10000;  % How many g samples to draw for the MC integration
factoryDefault.nDimsForMC = 1; % If there are this many dimensions, use MC integration
% 1 and 2D integration work fine with standard
% integral.m but 3d is slow and higher not
% possible.
factoryDefault.verbose = true; % Show messages
factoryDefault.useRal = false;  % Use the ral class (reprensent Reals As Logs) to improve numerical stability. 
%                       This is used only for anova and is necessary for
%                       problems with large number of observations. See rouderS.
factoryDefault.nrWorkers = 0;  % Set to zero to use serial execution (i.e. disable parfor loops in the code)   

%% Check what we already have, and update if needed

if ispref('bayesFactor') 
    machinePref =  getpref('bayesFactor','options');
else
    machinePref = struct; % Empty
end

%Assign any new factory defaullt to this machine
factoryUpdate =false;
fn =setdiff(fieldnames(factoryDefault),fieldnames(machinePref));
for i=1:numel(fn)
    machinePref.(fn{i}) = factoryDefault.(fn{i});
    factoryUpdate =true;
end

% Update the options based on user input
if nargin>0
    fn = fieldnames(opt);
    for i=1:numel(fn)
        machinePref.(fn{i}) = opt.(fn{i});
    end
end

if machinePref.useRal &&  machinePref.nDimsForMC >1
    warning('Using RAL only works for MC integration. Forcing MC.');
    machinePref.nDimsForMC =1;
end 



if ~ispref('bayesFactor') || nargin>0 || factoryUpdate
    % First run : setup pref
    % Or opt specified - save as pref for future use.
    setpref('bayesFactor','options',machinePref);
end
