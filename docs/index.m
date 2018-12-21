%% Bayes Factor 
% This is a Matlab class to run Bayes Factor analysis. 
% The code was written by Bart Krekelberg (bart@rutgers.edu) with some code
% taken from Sam Schwarzkopf's code and with inspriation from the R package
% by Richard Morey.
%
% The mathematical underpinning of these tests can be found in the
% following papers:
%
% # Rouder, J. N., Morey, R. D., Speckman, P. L. & Province, J. M. Default Bayes factors for ANOVA designs. J. Math. Psychol. 56, 356–374 (2012).
% # Kass, R. E. & Raftery, A. E. Bayes factors. J. Am. Stat. Soc. 90, 733–795 (1995).
% # Morey, R. D. & Wagenmakers, E. J. Simple relation between Bayesian order-restricted and point-null hypothesis tests. Stat. Probab. Lett. 92, 121–124 (2014).
%
% Currently the following statistical tests have been implemented
%
% * One sample t-test  (|ttest|)
% * Two sample t-test (|ttest2|)
% * N-Way Anova with fixed and random effects  (|linearMixedModel|)
% * Pearson Correlation  (|corr|)
% * Binomial Test  ( |binom| )
%
% For convenience, each of these tests also returns the results of a
% traditional, frequentist, test (i.e. p-values etc.)
%
%% Installation & Dependencies
% 
% The standard statistical analysis (and some utility functions) require the  
% Statistics and Machine Learning Matlab toolbox. 
% 
% All code is defined in a single class definition file (bayesFactor.m).
% Installation only requires adding the directory that contains this file to the
% Matlab search path. The |installBayesFactor| function does this for you.
%
%% Examples
% 