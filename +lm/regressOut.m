function T = regressOut(lme,pv)
% Regress out specific factors/effects specified by their full name as
% fixed, and all random effects (if random=true)
arguments
    lme (1,1)
    pv.random (1,1) = true  % Set to true to regress out all random effects
    pv.fixed (1,:) string  = ""  % Specify which effects to regress out
end

include = ~lme.ObservationInfo.Excluded;
T= lme.Variables(include,:);
% Regress out random effects
if pv.random
    Z = designMatrix(lme,'random');  % Random effects
    T.(lme.ResponseName)  = T.(lme.ResponseName) -Z(include,:)*lme.randomEffects;
end

% Regress out specified fixed effects
fe = lme.fixedEffects;
fe(~ismember(lme.Coefficients.Name,pv.fixed)) =0; 
if any(fe~=0)
    X =  designMatrix(lme,'fixed');
    T.(lme.ResponseName)  = T.(lme.ResponseName)-X(include,:)*fe; 
end
