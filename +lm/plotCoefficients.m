function [estimate,low,high] =plotCoefficients(m,fe,x,scale)
% Visualize the effects of a generalized linear model by pulling the betas
% back trough the inverse of the link function and displaying the effects
% as bars with error bars.
%
% INPUT
% lm - Linear model
% fe - Fixed effects to show. [All except the intercept]
% x - Horizontal axis. [1:nrFixedEffects].
% scale - Scale the effects to
%                      'Intercept' ; scale to the fixed intercept effect
%                       'RANDOMSTD' : scale to the standard deviation of
%                       the random effects.
%                       'none' - do not scale (default)
% OUTPUT
% estimate - Fxied effect estimates
% lower  - lower edge of CI
% upper - upper edge of CI
%
% BK - Feb 2020.

if nargin <4
    scale = 'none';
end

if nargin<2
    fe = '';
end
if isempty(fe)
    % Show all except the intercept
    allFe = m.CoefficientNames;%regexprep(lm.CoefficientNames,'_(?<level>[\w\d]+)\>','');
    out  = strcmpi(allFe,'(Intercept)');
    allFe(out) =[];
    [estimate,low,high] = lm.plotCoefficients(m,allFe,1:numel(allFe),scale); % recurse
    return
end

if ischar(fe)
    fe = {fe};
end

if isa(m,'GeneralizedLinearMixedModel')
    link  = m.Link.Inverse;
else
    link = @(x)(x);
end

coeffs = m.Coefficients;
nrFe = numel(fe);

if nargin <3 || isempty(x)
    x = 1:nrFe;
end

low = NaN(nrFe,1);
high = NaN(nrFe,1);
estimate = NaN(nrFe,1);
allFeNames = m.CoefficientNames; %regexprep(lm.CoefficientNames,'_(?<level>[\w\d]+)\>','');
for i=1:nrFe
    ix = find(ismember(allFeNames,fe{i}));
    low(i) = link(coeffs.Lower(ix));
    high(i) = link(coeffs.Upper(ix));
    estimate(i) =  link(coeffs.Estimate(ix));
end


switch upper(scale)
    case 'NONE'
        scale =1;
        units ='';
    case 'INTERCEPT'
        % Scale to the fixed effect intercept
        % How big is the effect relative to the grand mean.
        scale  = 0.01*m.Coefficients.Estimate(strcmpi(m.CoefficientNames,'(Intercept)'));
        units = '(% grand mean)';
    case 'RANDOMSTD'
        % Scale to the standard deviation of the random effects
        % "How big is the effect relative to the variation
        % across subjects?'
        scale = 0.01*std(m.randomEffects);
        units = '(% Random Effects Stdev)';
    otherwise
        error ('Unknown scale %s ', scale);
end

estimate=estimate/scale;
low = low/scale;
high= high/scale;

bar(x,estimate,'FaceColor','w','EdgeColor','k');
hold on
errorbar(x,estimate,estimate-low,high-estimate,'.','CapSize',0,'Color','k');
feNames =fe;% regexprep(fe,'_(?<level>[\w\d]+)\>','');
ylabel([m.ResponseName ' ' units])
set(gca,'XTick',x,'XTickLabel',strrep(feNames,'_','-'),'XtickLabelRotation',45);

end
