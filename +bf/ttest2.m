function [bf10,pValue] = ttest2(X,Y,varargin)
% Calculates Bayes Factor for a two-sample t-test.
% X - Sample 1
% Y - Sample 2 (not necessarily of the same size)
%
% EXAMPLES
% Provide samples X and Y
% [bf10,pValue] = bf.ttest2(X,Y,'tail','right')
% Or provide the results of the regular ttest2:
% [bf10,pValue]= bf.ttest2('T',T,'N',[10 20],'tail','both')
%
% Optional Parm/Value pairs:
% alpha - significance level for the frequentist test. [0.05]
% tail - 'both','right', or 'left' for two or one-tailed tests [both]
%               Note that  'right' means X>Y and 'left' is X<Y
% scale - Scale of the Cauchy prior on the effect size  [sqrt(2)/2]
% T     -  Instead of samples X,Y, provide the T and the number of samples in 
%           each of the two groups.
% N     - [Nx Ny] : samples in each of the groups. Can be a scalar for
%               equal samples.
% 
%
% OUTPUT
% bf10 - The Bayes Factor for the hypothesis that the means of the samples are different
% p     - p value of the frequentist hypothesis test
%
% Internally this code calls bf.ttest for the computation of Bayes Factors.
%
%
% BK - Nov 2018
% Nov 21 - changed interface to avoid confusion about df and N. Force
% specifying T and N. Compute df internally.


if isnumeric(X)
    parms = varargin;
else
    % Firs input is not numeric. 
    % % This must be a call with 'T' and 'N' specified
    parms = cat(2,{X,Y,},varargin);
    X=[];Y=[];
end

p=inputParser;
p.addParameter('alpha',0.05);
p.addParameter('tail','both',@(x) (ischar(x)&& ismember(upper(x),{'BOTH','RIGHT','LEFT'})));
p.addParameter('scale',sqrt(2)/2);
p.addParameter('T',[],@isnumeric);
p.addParameter('N',[],@isnumeric);
p.parse(parms{:});

if isempty(p.Results.T)
    % Calculate frequentist from the X and Y data
    tail = p.Results.tail;
    [~,~,~,stats] = ttest2(X,Y,'alpha',p.Results.alpha,'tail',tail);
    nX = numel(X);
    nY = numel(Y);
    T = stats.tstat;    
    N = [nX nY];
elseif isempty(p.Results.N)
    error('N must be specified when calling bf.ttest2 with a T');  
else
    T = p.Results.T;
    N = p.Results.N;
    if numel(N)==1
        %Equal samples
        N = [N N];
    end          
end

[bf10,pValue] = bf.ttest('T',T,'N',N,'scale',p.Results.scale,'tail',p.Results.tail);
end
