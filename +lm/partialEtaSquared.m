function [eta,lowerEta,upperEta] = partialEtaSquared(varargin)
% Returns the partial eta squared for each of the terms in the linear
% model, based on the ANOVA table.
% Assumes that each term in the model is a manipu;ated/controlled
% independent variable, and not a covariate
% INPUT
% 'lm' = Alinear model
% OR 
% 'F' = The F statistic
% 'df1' = Effect degrees of Freedom
% 'df2' - Error degrees of freedom.
% 'alpha' = Significance level for confidence intervals. [0.1] (i..e 90% CI
%           is default)
% 'nrSteps'  How many steps to use in the search for the bounds. This
%               determines the "resolution" of the bounds.  [1000]
%           There is a tradeoff with speed.
% OUTPUT
% v = The partial eta squared.
% lowerEta = The lower bound on the confidence interval
% upperEta = The upper bound on the confidence interval
% BK - Feb 2021.

p = inputParser;
p.addParameter('lm',[]);
p.addParameter('F',[]);
p.addParameter('df1',[]);
p.addParameter('df2',[]);
p.addParameter('alpha',0.1);
p.addParameter('nrSteps',1000);
p.parse(varargin{:});

if ~isempty(p.Results.lm)    
    F =p.Results.lm.anova.FStat;
    u = p.Results.lm.anova.DF1;
    v = p.Results.lm.anova.DF2;
else
    F = p.Results.F;
    u = p.Results.df1;
    v = p.Results.df2;
end

eta = F.*u./(F.*u+v);
FToDelta = @(f,U,V) (f.*(U./V).*(U+V+1)); 
if nargout>1    
    nrEtas = numel(eta);
    lowerEta = nan(size(eta));
    upperEta = nan(size(eta));
    delta = F.*(u./v).*(u+v+1);% Non centrality parameter delta
    for i=1:nrEtas        
        rangeOfDeltas =linspace(0,delta(i),p.Results.nrSteps); % Search up to delta.
        upperDelta = rangeOfDeltas(find(ncfinv((1-p.Results.alpha/2),u(i),v(i),rangeOfDeltas)>=F(i),1,'first'));
        lowerEta(i) = upperDelta./(upperDelta+u(i)+v(i)+1);
    end
end
if nargout >2
    for i=1:nrEtas
        % Estimate what the largest possible upperDelta will be ; that will
        % setup the search range. Not sure whether this always works. But
        % the assert will catch it if it doesn't.
        stopF = ncfinv(1-p.Results.alpha*0.0005,u(i),v(i),delta(i));
        stopDelta = stopF.*u(i)./v(i).*(u(i)+v(i)+1);
        rangeOfDeltas =linspace(delta(i),stopDelta,p.Results.nrSteps);
        upperDelta = rangeOfDeltas(find(ncfinv(p.Results.alpha/2,u(i),v(i),rangeOfDeltas)<=F(i),1,'last'));
        assert(upperDelta < max(rangeOfDeltas),'Search range for upper bound was too small. Change the heuristic...');            
        upperEta(i) = upperDelta./(upperDelta+u(i)+v(i)+1);
    end
end

