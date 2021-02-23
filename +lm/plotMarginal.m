function plotMarginal(m,base,effects,cumulative)
% Plot marginal effects for a (generalized) linear model, using a bar graph
% with error bars, and expressing all effects relative to a specified
% baseline.
% 
% INPUT 
% lm - linear model
% base - The base condition from which the marginal effects are
% computed.(I.e. this condition has effect zero and everything that is
% shown is relative to this). This should be a char.
% effects - Which (fixed) effects to show
% cumulative - Show effects as additional effects relative to the effects
% to the left in the graph.
%
% OUTPUT
% 
% BK - Feb 2020
if nargin <4
    cumulative = true;
end
coeffs = m.Coefficients;
nrFe = numel(effects);
delta = NaN(nrFe,1);
neg = NaN(nrFe,1);
pos= NaN(nrFe,1);

allFeNames = m.CoefficientNames; %regexprep(lm.CoefficientNames,'_(?<level>[\w\d]+)\>','');

if isa(m,'GeneralizedLinearMixedModel')
    link  = m.Link.Inverse;
else
    link = @(x)(x);
end
ix = ismember(allFeNames,base);
baseBeta = sum(coeffs.Estimate(ix));
baseEstimate=  link(baseBeta);

%% Compute
for i=1:numel(effects)
    ix = ismember(allFeNames,effects{i});
    
    thisUpper = link(baseBeta+coeffs.Upper(ix));
    thisEstimate = link(baseBeta+coeffs.Estimate(ix));
    thisLower =  link(baseBeta+coeffs.Lower(ix));
    
    delta(i) =  thisEstimate-baseEstimate;
    neg(i)   = thisEstimate-thisLower;
    pos(i) =  thisUpper-thisEstimate;    
    if cumulative
        baseBeta = baseBeta+coeffs.Estimate(ix);
        baseEstimate = thisEstimate;
    end
end
%% Visualize
x= 1:nrFe;
bar(x,delta,'FaceColor','w','EdgeColor','k');
hold on
errorbar(x,delta,neg,pos,'.','CapSize',0,'Color','k');
feNames =effects; regexprep(effects,'_(?<level>[\w\d]+)\>','');
yl = ylim;
[nudge,ix] = max(abs(yl));
nudge = 0.1*sign(yl(ix))*nudge;
text(-0.1,nudge,['Base = ' num2str(baseEstimate,2)],'HorizontalAlignment','Left');
ylabel(['\Delta  ' m.ResponseName ])
set(gca,'XTick',x,'XTickLabel',feNames);

end