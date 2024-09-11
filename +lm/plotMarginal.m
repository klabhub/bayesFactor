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
% to the left in the graph.[true]
%
% OUTPUT
% void
% BK - Feb 2020
arguments
    m (1,1) 
    base (1,:) char = '(Intercept)'
    effects (1,:) cell = {}
    cumulative (1,1) logical = true
end
coeffs = m.Coefficients;
if isempty(effects)
    effects= setdiff(m.CoefficientNames,{base}); % All except the base
end
nrFe = numel(effects);
delta = NaN(nrFe,1);
neg = NaN(nrFe,1);
pos= NaN(nrFe,1);

allFeNames = m.CoefficientNames; %regexprep(lm.CoefficientNames,'_(?<level>[\w\d]+)\>','');

if isa(m,'GeneralizedLinearMixedModel')
    unlink  =m.Link.Inverse;
else
    unlink = @(x)(x);
end
ix = ismember(allFeNames,base);
trueBaseBeta = sum(coeffs.Estimate(ix));
trueBaseEstimate=  unlink(trueBaseBeta);

%% Compute
for i=1:numel(effects)
    baseBeta = trueBaseBeta;
    if i>1 && cumulative
        % Add betas that have already been processed
        stay = ismember(allFeNames,effects(1:i-1));
        baseBeta = baseBeta + sum(coeffs.Estimate(stay));
    end
    baseEstimate = unlink(baseBeta);
     

     % Find the relevant effects 
    ix = ismember(allFeNames,effects{i});
    assert(any(ix),'%s does not exist in this model',effects{i});
    
    % Determine model estimates
    thisUpper = unlink(baseBeta+coeffs.Upper(ix));
    thisEstimate = unlink(baseBeta+coeffs.Estimate(ix));
    thisLower =  unlink(baseBeta+coeffs.Lower(ix));
    % Compare to base estimate
    delta(i) =  thisEstimate-baseEstimate;
    neg(i)   = thisEstimate-thisLower;
    pos(i)   =  thisUpper-thisEstimate;       
    
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
text(-0.1,nudge,sprintf('Base (%s) = %.2f',base,trueBaseEstimate),'HorizontalAlignment','Left','Interpreter','none');
ylabel(['\Delta  ' m.ResponseName ])
set(gca,'XTick',x,'XTickLabel',feNames);
ax=gca;
ax.XAxis.TickLabelInterpreter = 'none';
end