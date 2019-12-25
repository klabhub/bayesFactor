%% Example with continuous covariates (i.e. regression)

%% Liang et al.
% Compare the Liang et al formula (bf.bfFromR2) with the explciit
% integration implemented in the bf.anova function (so that it can handle
% continuous covariates as well as mixtures of continuous and categorical
% covariates).
nrSamples = 50;
% Create co-variates with different variance
nrRegressors = 2;
regressors = repmat([.1 3],[nrSamples 1]).*rand(nrSamples,nrRegressors);
noiseLevel = 1;    
nrIterations = 10; % Esimate N models
bfLiang = nan(nrIterations,1);
bfFull = nan(nrIterations,1);
for i=1:nrIterations
    beta = randn(1,nrRegressors);
    y = regressors*beta'+noiseLevel*randn([nrSamples 1]);
    [b,~,~,~,stats] = regress(y,[ones(nrSamples,1) regressors]);    
    R2 =stats(1);
    % Use the Liang et al formula directly
    bfLiang(i) = bf.bfFromR2(R2,nrSamples,nrRegressors);
    % Create a table to do integration in the bf.anova function
    r = num2cell(regressors,1);
    T = table(y,r{:},'VariableNames',{'y','x1','x2'});
    bfFull(i)=bf.anova(T,'y~x1+x2','continuousScale',1);
end
% Show the direct comparison. This should be a slope 1 line
figure(1);
clf;
subplot(1,2,1);
plot(bfLiang,bfFull,'.')
hold on
set(gca,'XScale','Log','YScale','Log')
plot(xlim,xlim,'k')
axis square;
xlabel 'BF Integrated'
ylabel 'BF Liang Formula'
subplot(1,2,2);
hist(bfFull./bfLiang)
xlabel 'BF Ratio (Anova/bfFromR2)'
ylabel '# Models'


%% Attitude
% Analyze the attitude data set  (See
% <https://richarddmorey.github.io/BayesFactor/#glm>)
load attitude
% Compare the rating~complaints model with a set of alternatives.
altModels = {'rating~complaints+learning','rating~complaints+learning+advance','rating~complaints+raises','rating~complaints+privileges','rating~complaints+advance'};
bfComplaints = bf.anova(attitude,'rating~complaints','alternativeModel',altModels,'continuousScale',sqrt(2)/4);
%% 
for i=1:numel(altModels)
    fprintf('%s:\t\t%3.2f\n',altModels{i},1./bfComplaints(i))
end

