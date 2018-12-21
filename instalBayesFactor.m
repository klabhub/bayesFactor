function installBayesFactor
% Install the bayesFactor class in Matlab.

lmm = which('linearmixedmodel'); 
if isempty(lmm)
    disp('Hmm. You don''t seem to have the the Statistics and Machine Learning Toolbox. This is not likely to work');
end

filename = mfilename('fullpath');
pth = fileparts(filename);
addpath(pth)
savepath;

end