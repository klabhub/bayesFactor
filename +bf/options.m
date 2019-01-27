function opt = options
%% Default parameters for the Monte Carlo integration and parallel execution
opt.minG = 0.0001;
opt.maxG = 10000;
opt.stepG = 0.05;
opt.nrSamples = 10000;
opt.nDimsForMC = 3; % If there are this many dimensions, use MC integration
% 1 and 2D integration work fine with standard
% integral.m but 3d is slow and higher not
% possible.
opt.nrWorkers = 4;  % Set to zero to use serial execution (i.e. disable parfor loops in the code)
end