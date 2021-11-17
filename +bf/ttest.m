function [bf10,pValue] = ttest(X,varargin)
%TTEST Bayes Factors for one-sample and paired t-tests.
% [bf10,p] = ttest(X,varargin)    - one sample
% [bf10,p] = ttest(X,M,varargin)   -one sample,non-zero mean
% [bf10,p] = ttest(X,Y,varargin)   -paired samples
%
% [bf10,p] = ttest('T',T,'N',10)   - calculate BF based
% on regular ttest output
%
% INPUT 
% X = single sample observations  (a column vector)
% Y = paired observations (column vector) or a scalar mean to compare the samples in X to.
%       [Defaults to 0]
%
% Optional Parm/Value pairs:
% tail - 'both','right', or 'left' for two or one-tailed tests [both]
% scale - Scale of the Cauchy prior on the effect size  [sqrt(2)/2]
%
% To calculated BF based on the outcome of a T-test, pass T and the number of
% samples as parm/value pairs:
% T - The T-value resulting from a standard T-Test output 
% N  - the number of samples used for the T-test
%
% OUTPUT
% bf10 - The Bayes Factor for the hypothesis that the mean is different
%           from zero. Using JZS priors. 
% p - p value of the frequentist hypothesis test
% CI    - Confidence interval for the true mean of X
%
% Based on: Rouder et al. J. Math. Psych. 2012
% 
% BK - Nov 2018
% Nov 21 - changed interface to avoid confusion about df and N. Force
% specifying T and N. Compute df internally.


if isnumeric(X)
    if mod(numel(varargin),2)==0
        % Only X specified
        Y = 0;
        parms = varargin;
    else
        % X and Y specified
        if numel(varargin)>1
            parms = varargin(2:end);
        else
            parms = {};
        end
        Y  =varargin{1};
    end
else
    %Neither X nor Y specified (must be a call with 'T' and 'N' specified
    parms = cat(2,X,varargin);
    X=[];Y=[];
end

p=inputParser;
p.addParameter('tail','both',@(x) (ischar(x)&& ismember(upper(x),{'BOTH','RIGHT','LEFT'})));
p.addParameter('scale',sqrt(2)/2);
p.addParameter('T',[],@isnumeric);
p.addParameter('N',[],@isnumeric);
p.parse(parms{:});

tail = p.Results.tail;

if isempty(p.Results.T)
    % Calculate frequentist from the X and Y data
    [~,pValue,~,stats] = ttest(X,Y,'tail',tail);
    T = stats.tstat;
    df = stats.df;
    N = numel(X);
elseif isempty(p.Results.N)
    error('N must be specified when calling bf.ttest with a T');  
else
    % User specified outcome of frequentist test (the builtin ttest);
    % Calculate BF from T and df.    
    T = p.Results.T;
    N = p.Results.N;
    if numel(N)==2
        % Call from a 2-sample T-tets (see bf.ttest2)
        % Adjust df and N
        df = sum(N)-2;
        N = prod(N)/sum(N);
    else
        df = N-1;
    end
    pValue = tcdf(T,df,'upper'); % Right tailed
        switch upper(tail)
            case 'BOTH'
                pValue = 2*(1-pValue );
            case 'LEFT'
                pValue   =1-pValue ;
            case 'RIGHT'
                % Ok as is
        end        
end

% Use the formula from Rouder et al.
% This is the formula in that paper; it does not use the
% scale numerator = (1+T.^2/(N-1)).^(-N/2);
% fun  = @(g) ( ((1+N.*g).^-0.5) .* (1+T.^2./((1+N.*g).*(N-1))).^(-N/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );
% Here we use the scale  (Checked against Morey's R package and
% http://pcl.missouri.edu/bayesfactor)
r = p.Results.scale;
numerator = (1+T.^2/df).^(-(df+1)/2);
fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

% Integrate over g
bf01 = numerator/integral(fun,0,inf);
% Return BF10
bf10 = 1./bf01;

switch (tail)
    case 'both'
        % Nothing to do
    case {'left','right'}
        % Adjust the BF using the p-value as an estimate for the posterior
        % (Morey & Wagenmakers, Stats and Prob Letts. 92 (2014):121-124.
        bf10 = 2*(1-pValue)*bf10;
end
end

