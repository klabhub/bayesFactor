function v = dummyVarCoding(m)
% Use as m.DummyVarCoding.
% 
% Extract dummy variable coding from the linear model. For some
% inexplicable reason this is a protected property of the linear models,
% sow e have to convert to struct first (and disable the associated warning
% temporarily).  Once Matlab makes DummyVarCoding public in their classreg implementation,
% we can remove this function and the calling code should keep working.
%
% INPUT
% m - linear mixed model
% OUTPUT
% v - string that identifes the dummy variable coding used in the model
%
% BK - Mar 2021


state = warning('query','MATLAB:structOnObject');
warning('off',   'MATLAB:structOnObject');
s = struct(m);
v = s.DummyVarCoding;
warning(state);