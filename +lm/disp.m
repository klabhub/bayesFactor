function out = disp(m,factors,showEffects,tol)
% Convenience disp function to show Anova results of a linear model in
% standard notation for easy copy and paste.
%
% INPUT
%  m =  a linear model
% factors = a cell array of factors whose stats are to be shown. Defaults
% to all but the Intercept.
% showEffects = Show effect sizes.
%                       PARTIALETA - an standardized effect size , may be
%                       of limited use for LMM
%                      'Intercept' ; scale to the fixed intercept effect
%                       'RANDOMSTD' : scale to the standard deviation of
%                       the random effects.
%                       'RAW' - do not scale
% tol = Tolerance for partial eta squared confidence intervals
% OUTPUT
% out = the string that is also written to the command line .
%
% BK - Feb 2020
if nargin<2 || isempty(factors)
    factors =m.anova.Term(2:end);
end
if nargin <3
    showEffects = 'raw';
end
if islogical(showEffects)
    showEffects = 'raw';
end
if nargin <4
    tol =0.001;
end
if nargout>0
    out = '';
end

if ischar(factors)
    if strcmpi(factors,'*')
        factors =m.anova.Term; % All including the intercept
    else
        factors = {factors};
    end
end

fprintf('%s\n',m.Formula.char);
if contains(showEffects,'PARTIALETA','IgnoreCase',true)
    [partialEta,partialEtaLB,partialEtaUB] = lm.partialEtaSquared(m,'tol',tol);
    eta = ['partial ' char(hex2dec('03B7')) char(178)];
    hasEta = true;
else
    hasEta =false;
end

for f=1:numel(factors)
    factor = factors{f};
    stay = strcmpi(m.anova.Term,factor);
    
    fmt = '\t %s: F(%d,%d)=%3.1f, p=%3.2g,';
    vars = {factor,m.anova.DF1(stay,1),m.anova.DF2(stay,1),m.anova.FStat(stay,1),m.anova.pValue(stay,1)};
    
    if hasEta
        fmt = [fmt eta '=%3.3g CI: [%4.4g, %4.4g]'];     %#ok<AGROW>
        vars = cat(2,vars,{partialEta(stay),partialEtaLB(stay),partialEtaUB(stay)});
    end
    
    if m.anova.pValue(stay,1) <0.05
        style = 2; % Error output stream ; red
    else
        style =1; % stdout;
    end

    [scale,units] = lm.scaleFactor(m,showEffects);
    
    fe = m.fixedEffects;
    elms = strsplit(factor,':');
    nrElms =numel(elms);
    for e=1:nrElms
        if strcmpi(elms{e},'(Intercept)')
            % not in variable info (and not categorical, so continue
            elms{e} = '\(Intercept\)';
            if strcmpi(showEffects,'Intercept')
                % Revert to RAW scaling for the intercept itself.
                scale = 1;
                units = '';
            end
        elseif m.VariableInfo{elms{e},'IsCategorical'}
            elms{e} = [elms{e} '_[^:]+'];
        end
        if nrElms>1 && e <nrElms
            elms{e} = [elms{e} ':'];
        end
    end
    expression  = ['^' strcat(elms{:}) '$'] ;
    stayFe = ~cellfun(@isempty, regexp(m.CoefficientNames,expression));
    
    % Find the corresponding linear effect (or the largest one for
    % a multilevel categorical factor).
    
    fe = fe(stayFe)/scale;
    [fe,ix] = max(fe);
    low = m.Coefficients.Lower(stayFe);
    low = low(ix)/scale;
    up =m.Coefficients.Upper(stayFe);
    up = up(ix)/scale;
    vars = cat(2,vars,{fe,low,up});
    if sum(stayFe)>1
        fmt = cat(2, fmt ,[' (Max Effect: %3.3g ' units ' CI [%4.4g, %4.4g])']);
    else
        fmt = cat(2, fmt ,[' (Effect: %3.3g ' units '  CI [%4.4g, %4.4g])']);
    end
    
    fprintf(style,[fmt '\n'],vars{:});
    if nargout>0
        out = char(out,sprintf([fmt '\n'],vars{:}));
    end
end
end
