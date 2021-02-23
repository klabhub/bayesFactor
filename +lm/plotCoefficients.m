function [estimate,lower,upper] =plotCoefficients(m,fe,x)
% Visualize the effects of a generalized linear model by pulling the betas
% back trough the inverse of the link function and displaying the effects
% as bars with error bars. 
%
% INPUT
% lm - Linear model
% fe - Fixed effects to show. [All except the intercept]
% x - Horizontal axis. [1:nrFixedEffects].
% OUTPUT
% estimate - Fxied effect estimates 
% lower  - lower edge of CI
% upper - upper edge of CI
% 
% BK - Feb 2020.

if nargin<2
    fe = '';
end
if isempty(fe)
    % Show all except the intercept
    allFe = m.CoefficientNames;%regexprep(lm.CoefficientNames,'_(?<level>[\w\d]+)\>','');
    out  = strcmpi(allFe,'(Intercept)');
    allFe(out) =[];
    [estimate,lower,upper] = m.plotCoefficients(m,allFe,1:numel(allFe));
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

if nargin <3
    x = 1:nrFe;
end

lower = NaN(nrFe,1);
upper = NaN(nrFe,1);
estimate = NaN(nrFe,1);
allFeNames = m.CoefficientNames; %regexprep(lm.CoefficientNames,'_(?<level>[\w\d]+)\>','');
for i=1:nrFe
    ix = find(ismember(allFeNames,fe{i}));
    lower(i) = link(coeffs.Lower(ix));
    upper(i) = link(coeffs.Upper(ix));
    estimate(i) =  link(coeffs.Estimate(ix));
    
end
bar(x,estimate,'FaceColor','w','EdgeColor','k');
hold on
errorbar(x,estimate,estimate-lower,upper-estimate,'.','CapSize',0,'Color','k');
feNames =fe;% regexprep(fe,'_(?<level>[\w\d]+)\>','');
ylabel(m.ResponseName)
set(gca,'XTick',x,'XTickLabel',feNames);

end
