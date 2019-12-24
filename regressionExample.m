%% Example with continuous covariates (i.e. regression)

% Compare the Liang et al formula (bf.bfFromR2) with the explciit
% integration implemented in the anova function (so that it can handle
% continuous covariates as well as mixtures of continuous and categorical
% covariates).
nrSamples = 50;
nrRegressors =3;
regressors = rand(nrSamples,nrRegressors);
noiseLevel = 2.0;   
nrIterations = 10;
bfLiang = nan(nrIterations,1);
bfFull = nan(nrIterations,1);
for i=1:nrIterations
    beta = randn(1,nrRegressors);
    y = regressors*beta'+noiseLevel*randn([nrSamples 1]);
    [b,~,~,~,stats] = regress(y,[ones(nrSamples,1) regressors]);    
    R2 =stats(1);
    % Use the Liang et al formula directly
    bfLiang(i) = bf.bfFromR2(R2,nrSamples,nrRegressors);
    % Create a table to do integration in the anova function
    r = num2cell(regressors,1);
    T = table(y,r{:},'VariableNames',{'y','x1','x2','x3'});
    bfFull(i)=bf.anova(T,'y~x1+x2+x3');
end
%%
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
hist(log10(bfFull)-log10(bfLiang))
xlabel 'Difference (log10)'
ylabel '# Models'