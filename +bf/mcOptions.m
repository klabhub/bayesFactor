function opt = mcOptions
  %% Default parameters for the Monte Carlo integration.
        opt.minG = 0.0001;
        opt.maxG = 10000;
        opt.stepG = 0.05;
        opt.nrSamples = 10000;
        opt.nDimsForMC = 3; % If there are this many dimensions, use MC integration
        % 1 and 2D integration work fine with standard
        % integral.m but 3d is slow and higher not
        % possible.
end