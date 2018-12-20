function [Qa,Xa] = fixedEffectConstraint(X)

if isscalar(X)
    % Thisis the number of effects
    nrEffects = X;
else
    %Design matrix was passed
    nrEffects = size(X,2);
end

Sigmaa =eye(nrEffects)- ones([nrEffects nrEffects])/nrEffects;
[eigenVecs,ev]= eig(Sigmaa','vector');
[~,ix] = sort(ev,'desc');
Qa = eigenVecs(:,ix(1:end-1));
%Iaminus1 = eye(nrEffects-1);
%Sigmaa = Qa*Iaminus1*Qa';
if isscalar(X)
    Xa =[];
else
    % Transform design matrix
    Xa = X*Qa;
end

