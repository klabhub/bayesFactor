function [eta,lowerEta,upperEta] = partialEtaSquared(m,varargin)
% Returns the partial eta squared for each of the terms in the linear
% model, based on the ANOVA table.
%
% Assumes that each term in the model is a manipu;ated/controlled
% independent variable, and not a covariate
% 
% INPUT
% 'lm' = A linear model
% OR 
% 'F' = The F statistic
% 'df1' = Effect degrees of Freedom
% 'df2' - Error degrees of freedom.
% 'alpha' = Significance level for confidence intervals. [0.1] (i..e 90% CI
%           is default)
% 'tol'  =  Confidence bounds are determined with this tolerance on the bounds. 
%           The default is 0.005 which is in the same units as eta, so corresponds to 
%           a tolerance of 0.5% variance explained. Setting this smaller slows the
%           computation down.
% OUTPUT
% v = The partial eta squared.
% lowerEta = The lower bound on the confidence interval
% upperEta = The upper bound on the confidence interval
% BK - Feb 2021.

p = inputParser;
p.addOptional('m',[],@(x)isa(x,'classreg.regr.LinearLikeMixedModel'));
p.addParameter('F',[]);
p.addParameter('df1',[]);
p.addParameter('df2',[]);
p.addParameter('alpha',0.1);
p.addParameter('tol',0.005);
p.parse(m,varargin{:});

if ~isempty(p.Results.m)    
    F =p.Results.m.anova.FStat;
    u = p.Results.m.anova.DF1;
    v = p.Results.m.anova.DF2;
else
    F = p.Results.F;
    u = p.Results.df1;
    v = p.Results.df2;
end


eta = F.*u./(F.*u+v);


if nargout>1    
    nrEtas = numel(eta);
    lowerEta = nan(size(eta));
    upperEta = nan(size(eta));    
    for i=1:nrEtas        
        etaRange = 0:p.Results.tol:eta(i);  % This is the range of etas that we will test for LB: up to eta in steps of tol
        if numel(etaRange)==1 % 0 is within tol of eta(i)
            lowerEta(i) = 0;
        else           
            rangeOfDeltas = etaRange.*(u(i)+v(i)+1)./(1-etaRange);        % Convert to delta (noncentrality parameter)        
            upperDelta = rangeOfDeltas(find(ncfinv((1-p.Results.alpha/2),u(i),v(i),rangeOfDeltas)>=F(i),1,'first')); % Find alpha%        
            lowerEta(i) = upperDelta./(upperDelta+u(i)+v(i)+1); % Convert back to eta
        end
    end
end
if nargout >2
    for i=1:nrEtas
        etaRange = eta(i):p.Results.tol:(1-p.Results.tol); % Search from eta and upward
        if numel(etaRange) ==1 % eta is within tol of 1
            upperEta(i) =  1;
        else
            rangeOfDeltas = etaRange.*(u(i)+v(i)+1)./(1-etaRange);
            upperDelta = rangeOfDeltas(find(ncfinv(p.Results.alpha/2,u(i),v(i),rangeOfDeltas)<=F(i),1,'last'));
            if isempty(upperDelta)
                %For very small F this happens: weird, even eta(i) is not small enough... so the upper
                %limit is below eta(i)...?
                upperEta(i) = NaN;
            else
                assert(upperDelta < max(rangeOfDeltas),'Search range for upper bound was too small. Change the heuristic...');            
                upperEta(i) = upperDelta./(upperDelta+u(i)+v(i)+1);
            end
        end
    end
end


