classdef bayesFactor < handle
    % A class to perform Bayes Factor statistical analysis to quantify
    % evidnce in favor or against a hypothesis.
    % For background see:
    %
    % The mathemetical underpinning for these routines is provided in
    %
    %
    % Implementation notes:
    %  The class currently does not store data, only some default parameers
    %  that are used in different contexts.
    % BK - 2018
    
    properties (SetAccess = public, GetAccess =public)
        %% Default parameters for the Monte Carlo integration.
        minG = 0.0001;
        maxG = 10000;
        stepG = 0.05;
        nrSamples = 10000;
        nDimsForMC = 3; % If there are this many dimensions, use MC integration
        % 1 and 2D integration work fine with standard
        % integral.m but 3d is slow and higher not
        % possible.
    end
    
    %% Public acces functions - the main interface
    methods (Access=public)
        
        function o = bayesFactor
            % Constructor. Nothing to do.
            
        end
        
        function [bf10,model,aov] = linearMixedModel(o,tbl,formula,varargin)
            % Analyze table data with a linear mixed effects model.
            %
            % INPUT
            % tbl = The table with data.
            % formula = The linear (mixed) model using Wilcoxon notation.
            %           See LinearModel and LinearMixedModel in the stats
            %           toolbox for definition
            % Parm/Value pairs:
            % 'sharedPriors' - Which columns (i.e. factors) in the table should share a
            %                   their prior on their effect size. The default is that
            %                   all levels within a factor share a prior.
            %                   To share priors across factors, use
            %                      {'{'a','b'},{'c','d'}}   -> share priors for a and b
            %                       and, separately, for c and d.
            %                      {{'a','b','c','d'}} ->  share priors for all factors
            %                      (this is the 'single g' approach in Rouder.
            %                      Shortcuts:
            %                       'within' - share within a fixed effect factor, not across
            %                       'singleG' - share across all fixed effects.
            % 'interactions'    - Include interactions ('all') or not ('none'). Default
            %                       is  'none'.  Or specify a set {'a:b','c:d'}
            %
            % OUTPUT
            % bf10 - The Bayes Factor comparing the model to the model with intercept only.
            %       To compute BF for more refined hypotheses you compute
            %       a BF for the full model, and a restricted model and
            %       then take the ratio. See rouderFigures for examples.
            % model - A linear model or linear mixed model from the
            %           statistics toolbox
            % aov    - results of an ANOVA
            %
            % BK -2018
            
            p=inputParser;
            p.addParameter('interactions','none',@(x) (ischar(x) || iscell(x)));
            p.addParameter('sharedPriors','within',@(x) ischar(x) || (iscell(x) && iscell(x{1}))); % Cell containing cells with factors(columns) that share a prior.
            p.addParameter('treatAsRandom',{});
            p.parse(varargin{:});
            
            formula = classreg.regr.LinearMixedFormula(formula);
            if ~isempty(formula.RELinearFormula)
                 error('Sorry, grouping has not been implemented yet.)');
            else
                feFormula = formula.FELinearFormula;
                isMain= ~cellfun(@(x) (contains(x,'(')|| contains(x,':')),feFormula.TermNames);
                mainEffects  =feFormula.TermNames(isMain)';
                isInteraction= cellfun(@(x) (contains(x,':')),feFormula.TermNames);
                interactions = feFormula.TermNames(isInteraction)';
                response = formula.ResponseName;
            end
            nrMainEffects =numel(mainEffects);
            nrInteractions = numel(interactions);
            allTerms = cat(2,mainEffects,interactions);
            nrAllTerms = nrMainEffects+nrInteractions;
            
            %% Setup sharing of priors as requested
            if ischar(p.Results.sharedPriors)
                switch upper(p.Results.sharedPriors)
                    case 'WITHIN'
                        % Share priors for each level of each factor, but not across factors
                        sharedPriors = cell(1,nrAllTerms);
                        [sharedPriors{:}] = deal(allTerms{:});
                    case 'SINGLEG'
                        sharedPriors = {allTerms};
                    case {'NONE',''}
                        sharedPriors ={};
                end
            else
                sharedPriors = p.Results.sharedPriors;
            end
            sharedPriorIx = cell(1,numel(sharedPriors));
            [sharedPriorIx{:}] = deal([]);
            
            %% Construct the design matrix for the main effects
            designMatrix = cell(1,nrMainEffects +nrInteractions);
            for i=1:nrMainEffects
                thisX = classreg.regr.modelutils.designmatrix(tbl,'model','linear','intercept',false,'DummyVarCoding','full','PredictorVars',mainEffects{i},'responseVar','');
                if ismember(mainEffects{i},p.Results.treatAsRandom)
                    % Random effect - Keep full dummy X
                else
                    % Sum-to-zero contrasts that equates marginal priors across levels.
                    thisX = o.zeroSumConstraint(thisX);
                end
                designMatrix{i} =thisX; % Store for later use
            end
            
            %% Construct the interaction parts of the design matrix for all requested interactions
            for i=1:nrInteractions
                aName =extractBefore(interactions{i},':');
                bName = extractAfter(interactions{i},':');
                thisA = classreg.regr.modelutils.designmatrix(tbl,'model','linear','intercept',false,'DummyVarCoding','full','PredictorVars',aName,'responseVar','');
                thisB = classreg.regr.modelutils.designmatrix(tbl,'model','linear','intercept',false,'DummyVarCoding','full','PredictorVars',bName,'responseVar','');
                if ~ismember(aName,p.Results.treatAsRandom)
                    thisA = o.zeroSumConstraint(thisA);
                end
                if ~ismember(bName,p.Results.treatAsRandom)
                    thisB = o.zeroSumConstraint(thisB);
                end
                thisX = o.interaction(thisA,thisB);
                %thisX = thisX -mean(thisX,2);
                designMatrix{nrMainEffects+i} = thisX;
            end
            
            
            %% Assign groups of effects to use the same prior on effect size.
            soFar  =0;
            if isempty(sharedPriors)
                sharedPriorIx = {};
            else
                for i=1:numel(allTerms)
                    match = cellfun(@any,cellfun(@(x)(strcmp(allTerms{i},x)),sharedPriors,'UniformOutput',false));
                    if ~any(match)
                        error(['Shared priors not defined for ' allTerms{i}]);
                    end
                    nrInThisTerm  = size(designMatrix{i},2);
                    sharedPriorIx{match} = cat(2,sharedPriorIx{match},soFar+(1:nrInThisTerm));
                    soFar = soFar+nrInThisTerm;
                end
            end
            %% Call the anova function for the actual analysis
            bf10 = o.nWayAnova(tbl.(response),[designMatrix{:}],'sharedPriors',sharedPriorIx);
            
            
            if nargout>1
                % Traditional
            end
        end
        
        
        
    end
    
    %% Internal computations.
    methods (Access = protected)
        
        function v = mcIntegral(o,fun,prior,nrDims)
            % Monte Carlo integration
            %
            % INPUT
            % fun -  The function to integrate. This should be specified as a
            %       function_handle that takes a single input (g)
            % prior - the prior distribution of the g's. A function_handle.
            %
            % nrDims - The number of dimensions to integrate over. [1]
            % options - A struct with options.  []
            % OUTPUT
            % v -  The value of the integral. (Typically the BF10).
            
            
            %% Setup the PDF to do importance sampling
            gRange =  o.minG:o.stepG:o.maxG;
            pdf = prior(gRange);
            pdf = pdf./sum(pdf);
            % Draw samples weighted by this prior.
            g =nan(nrDims,o.nrSamples);
            for d=1:nrDims
                g(d,:) = randsample(gRange,o.nrSamples,true,pdf);
            end
            %% Evaluate the function at these g values
            bf10Samples = fun(g);
            pg = prod(prior(g),1);  % Probability of each g combination
            v = mean(bf10Samples./pg); % Expectation value- = integral.
        end
        
        
        function bf10 = nWayAnova(o,y,X,varargin)
            % ANOVA BF
            % y = data values
            % X = design matrix for ANOVA (indicator vars)  no constant term
            %
            % Parm/Value pairs
            % 'sharedPriors'  - Cell array of vectors indicating which effects (columns
            % of X) share the same prior. [{1:nrEffects}]: all effects share the same prior.
            %
            % BK 2018
            nrEffects = size(X,2);
            
            p =inputParser;
            p.addParameter('sharedPriors',{},@iscell); % Which effects share a prior? A cell array with indices corresponding to columns of X
            p.parse(varargin{:});
            
            if isempty(p.Results.sharedPriors)
                sharedPriors = {1:nrEffects};
            else
                sharedPriors = p.Results.sharedPriors;
            end
            
            prior = @(g)(bayesFactor.scaledInverseChiPdf(g,1,1));
            integrand = @(varargin) (bayesFactor.rouderS(cat(1,varargin{:}),y,X,sharedPriors).*prod(prior(cat(1,varargin{:})),1));
            nrDims = numel(sharedPriors);
            if nrDims>= o.nDimsForMC
                % Use MC Sampling to calculate the integral
                bf10 = o.mcIntegral(integrand,prior,nrDims);
            else
                switch (nrDims)
                    case 1
                        bf10 = integral(integrand,0,Inf);
                    case 2
                        bf10 = integral2(integrand,0,Inf,0,Inf);
                    case 3
                        bf10 = integral3(integrand,0,Inf,0,Inf,0,Inf);
                end
            end
        end
        
        
    end
    
    %% Helper functions
    methods (Static, Hidden)
        function G = gMatrix(grouping,g)
            % Generate a matrix in which each row corresponds to an effect, each column a value
            % that will be integrated over. The grouping cell array determines which of
            % the effects share a prior (i.e. levles of the same factor) and which have
            % their own.
            % grouping   - Cell array with vectors that contain effect indices (i.e.
            % columns of the design%matrix) that share a prior.
            % g  - The values for each of the priors. Each row is an independent prior,
            %               each column is a sample
            %
            % EXAMPLE
            % gMatrix({1 2],[3 4]},[0 0.1 0.2 0.3; 0.6 0.7 0.8 0.9])
            % g = [ 0 0.1 0.2 0.3;
            %       0 0.1 0.2 0.3;
            %       0.6 0.7 0.8 0.9;
            %       0.6 0.7 0.8 0.9]
            
            
            nrEffects = sum(cellfun(@numel,grouping));
            assert(nrEffects>0,'The number of groups must be at least one');
            nrValues = size(g,2);
            G = nan(nrEffects,nrValues);
            for i=1:numel(grouping)
                G(grouping{i},:) = repmat(g(i,:),[numel(grouping{i}) 1]);
            end
        end
        
        function value= rouderS(g,y,X,grouping)
            % The S(g) function of Eq 9 in Rouder et al.
            % g = Matrix of g values, each row is an effect, each column is a value
            % that we're integrating over.
            % y = Data values
            % X = design matrix, without a constant term, with indicator variables only
            
            g = bayesFactor.gMatrix(grouping,g);
            
            nrObservations = size(X,1);
            one = ones(nrObservations,1);
            P0 = 1./nrObservations*(one*one');
            yTilde = (eye(nrObservations)-P0)*y;
            XTilde = (eye(nrObservations)-P0)*X;
            nrPriorValues=size(g,2);
            value= nan(1,nrPriorValues);
            for i=1:nrPriorValues
                if all(g(:,i)==0)
                    value(i)=0;
                else
                    G = diag(g(:,i));
                    invG = diag(1./g(:,i));
                    Vg = XTilde'*XTilde + invG;
                    yBar = one'*y/nrObservations;
                    preFactor= 1./(sqrt(det(G))*sqrt(det(Vg)));
                    numerator =    y'*y-nrObservations*yBar^2;
                    denominator = ((yTilde'*yTilde) -yTilde'*XTilde*(Vg\XTilde'*yTilde));
                    value(i)= preFactor*(numerator/denominator).^((nrObservations-1)/2);
                end
            end
        end
        
        
        function [Xa,Qa] = zeroSumConstraint(X)
            % Impose a zero-sum constraint on a dummy coded predictor matrix X
            % as in Rouder et al. 2012 . The goal is to make the model estimable
            % while equating marginal priors across levels.
            %
            % By default this is applied to all fixed effects.
            %
            % INPUT
            % X = Dummy coded predictor Matrix [nrObservations nrLevels].
            %       By passing a scalar, the function computes only the projection
            %       matrix (Qa) for a predictor matrix with that many effects.
            % OUTPUT
            % Xa = Matrix with zero-sum constraint [nrObservations nrLevels-1]
            % Qa -  projection matrix (Xa = X*Qa).
            %
            % BK - 2018
            
            if isscalar(X)
                % This is the number of effects
                nrEffects = X;
            else
                %Design matrix was passed
                nrEffects = size(X,2);
            end
            
            %% Follow  Rouder et al. 2012
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
            
        end
        
        function interactions = allInteractions(factorNames)
            % Given a cell array of factor names, create all pairwise combinations
            % of factors to represent interactions
            % INPUT
            % factorNames  - cell array of names
            % OUTPUT
            % interactionNames - cell array of interaction names
            %
            % allInteractions({'a','b'}) -> {'a:b'}
            %
            % BK  -Nov 2018
            
            cntr=0;
            nrFactors = numel(factorNames);
            nrInteractions = nrFactors*(nrFactors-1)-1;
            interactions =cell(1,nrInteractions);
            for i=1:nrFactors
                for j=(i+1):nrFactors
                    cntr= cntr+1;
                    interactions{cntr} = [factorNames{i} ':' factorNames{j}];
                end
            end
        end
        
        function X = interaction(Xa,Xb)
            % Create all interaction terms from two dummy coded design matrices.
            % (See Box II in Rouder et al. 2012)
            %
            % INPUT
            % Xa, Xb = [nrObservations nrA] and [nrObservations nrB]
            %           design matrices with matching number of observation (rows)
            % OUTPUT
            % X = Dummy coded design matrix [nrObservations nrA*nrB]
            %
            nA = size(Xa,2);
            nB = size(Xb,2);
            Xa = repmat(Xa,[1 nB]);
            Xb = repmat(Xb,[1 nA]);
            ix = (repmat(1:nA:nA*nB,[nA 1]) + repmat((0:nA-1)',[1 nB]))';
            Xa = Xa(:,ix(:));
            X= Xa.*Xb;
            
        end
        
        
        function y = inverseGammaPdf(x,alpha,beta)
            % The inverse Gamma PDF.
            % INPUT
            % x  (>0)
            % alpha - shape parameter
            % beta  - scale parameter
            %
            % BK - 2018
            %assert(all(x>0),'The inverse gamma PDF is only defined for x>0')
            z = x<0;
            y = zeros(size(z));
            y(~z) = (beta.^alpha)/gamma(alpha)*(1./x(~z)).^(alpha+1).*exp(-beta./x(~z));
        end
        
        function y = scaledInverseChiPdf(x,df,scale)
            % The Scaled Inverse Chi-Squared Distribution.
            % INPUT
            % x = the parameter value
            % df = degrees of freedom
            % scale = scale (tau squared).
            % OUTPUT
            % y = The probaility density
            %
            % BK - 2018
            %assert(all(x>0),'The scaled inverse Chi-squared PDF is only defined for x>0')
            z = x<0;
            y = zeros(size(z));
            if nargin <3
                scale =1; % Default to scaled inverse Chi-squared.
            end
            y(~z) = bayesFactor.inverseGammaPdf(x(~z),df/2,df*scale/2);
        end
        
    end
    
    %% Public User interface
    methods (Static, Hidden=false)
        function [bf10,p,CI,stats] = ttest2(X,Y,varargin)
            % Calculates Bayes Factor for a two-sample t-test.
            % X - Sample 1
            % Y - Sample 2 (not necessarily of the same size)
            %
            % Optional Parm/Value pairs:
            % alpha - significance level for the frequentist test. [0.05]
            % tail - 'both','right', or 'left' for two or one-tailed tests [both]
            % scale - Scale of the Cauchy prior on the effect size  [sqrt(2)/2]
            % stats - A struct containing .tstat  , .df , .pvalue .tail and .N - This allows one to
            %               calculate BF10 directly from the results of a standard ttest2 output.
            %           Note, however, that .N should be adjusted to nX*nY/(nX+nY) ,
            %           and                 .df = nx+ny-2
            %           If you call this fuction with data (X, Y) this adjustment is
            %           done automatically.
            %
            % OUTPUT
            % bf10 - The Bayes Factor for the hypothesis that the means of the samples are different
            % p     - p value of the frequentist hypothesis test
            % CI    - Confidence interval for the true mean of X
            % stats - Structure with .tstat, .df,resulting from the traditional test.
            %
            % Internally this code calls bf.ttest for the computation of Bayes Factors.
            %
            %
            % BK - Nov 2018
            
            if isnumeric(X)
                parms = varargin;
            else
                %Neither X nor Y specified (this must be a call with 'stats' specified
                parms = cat(2,{X,Y,},varargin);
                X=[];Y=[];
            end
            
            p=inputParser;
            p.addParameter('alpha',0.05);
            p.addParameter('tail','both',@(x) (ischar(x)&& ismember(upper(x),{'BOTH','RIGHT','LEFT'})));
            p.addParameter('scale',sqrt(2)/2);
            p.addParameter('stats',[],@isstruct);
            p.parse(parms{:});
            
            if isempty(p.Results.stats)
                % Calculate frequentist from the X and Y data
                tail = p.Results.tail;
                [~,p,CI,stats] = ttest2(X,Y,'alpha',p.Results.alpha,'tail',tail);
                nX = numel(X);
                nY = numel(Y);
                statsForBf = stats;
                statsForBf.p = p;
                statsForBf.N = nX*nY/(nX+nY);
                statsForBf.df = nX+nY-2;
                statsForBf.tail = tail;
            else
                % User specified outcome of frequentist test (the builtin ttest), calculate BF from T and
                % df.
                statsForBf = p.Results.stats;
            end
            
            bf10 = bayesFactor.ttest('stats',statsForBf);
        end
        
        
        
        function [bf10,pValue,CI,stats] = ttest(X,varargin)
            % function [bf10,p,CI,stats] = ttest(X,Y,varargin)  - paired
            % function [bf10,p,CI,stats] = ttest(X,varargin)    - one sample
            % function [bf10,p,CI,stats] = ttest(X,M,varargin)   -one sample,non-zero mean
            %
            % Calculates Bayes Factor for a one-sample or paired T-test.
            %
            % X = single sample observations  (a column vector)
            % Y = paired observations (column vector) or a scalar mean to compare the samples in X to.
            %       {Defaults to 0]
            %
            % Optional Parm/Value pairs:
            % alpha - significance level for the frequentist test. [0.05]
            % tail - 'both','right', or 'left' for two or one-tailed tests [both]
            % scale - Scale of the Cauchy prior on the effect size  [sqrt(2)/2]
            % stats - A struct containing .tstat  , .df , .pvalue .tail and .N - This allows one to
            %               calculate BF10 directly from the results of a standard T-Test output.
            %
            % OUTPUT
            % bf10 - The Bayes Factor for the hypothesis that the mean is different
            %           from zero
            % p - p value of the frequentist hypothesis test
            % CI    - Confidence interval for the true mean of X
            % stats - Structure with .tstat, .df,
            %
            % BK - Nov 2018
            
            if isnumeric(X)
                if iseven(numel(varargin))
                    % Only X specified
                    Y = 0;
                    parms = varargin;
                else
                    % X and Y specified
                    parms = varargin{2:end};
                    Y  =varargin{1};
                end
            else
                %Neither X nor Y specified (must be a call with 'stats' specified
                parms = cat(2,X,varargin);
                X=[];Y=[];
            end
            p=inputParser;
            p.addParameter('alpha',0.05);
            p.addParameter('tail','both',@(x) (ischar(x)&& ismember(upper(x),{'BOTH','RIGHT','LEFT'})));
            p.addParameter('scale',sqrt(2)/2);
            p.addParameter('stats',[],@isstruct);
            p.parse(parms{:});
            
            
            if isempty(p.Results.stats)
                % Calculate frequentist from the X and Y data
                tail = p.Results.tail;
                [~,pValue,CI,stats] = ttest(X,Y,'alpha',p.Results.alpha,'tail',tail);
                T = stats.tstat;
                df = stats.df;
                N = numel(X);
            else
                % User specified outcome of frequentist test (the builtin ttest), calculate BF from T and
                % df.
                T = p.Results.stats.tstat;
                df = p.Results.stats.df;
                pValue = p.Results.stats.p;
                tail  = p.Results.stats.tail;
                N = p.Results.stats.N;
                CI = [NaN NaN];
            end
            
            % Use the formula from Rouder 2009
            r = p.Results.scale;
            numerator = (1+T.^2/df).^(-(df+1)/2);
            fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );
            % Integrate over g
            bf01 = numerator/integral(fun,0,inf);
            % Return BF10
            bf10 = 1./bf01;
            
            switch (tail)
                case 'both'
                    % Nothing to do
                case {'left','right'}
                    % Adjust the BF using hte p-value as an estimate for the posterior
                    % (Morey & Wagenmakers, Stats and Prob Letts. 92 (2014):121-124.
                    bf10 = 2*(1-p)*bf10;
            end
        end
        
        function [bf10,p,pHat] = binotest(k,n,p)
            % Bayes factor for binomial test with k successes, n trials and base probability p.
            % INPUT
            %  k - number of successes
            %  n - number of draws
            %  p - true binomial probabiliy
            % OUTPUT
            % bf - Bayes Factor representing the evidence that this n/k
            % could result from random draws with p (BF>1) or not (BF<1)
            % p - p-value of a traditional test
            % pHat - esttimae of the binomial probablity
            
            % Code from Sam Schwarzkopf
            F = @(q,k,n,p) nchoosek(n,k) .* q.^k .* (1-q).^(n-k);
            bf01 = (nchoosek(n,k) .* p.^k .* (1-p).^(n-k)) / integral(@(q) F(q,k,n,p),0,1);
            bf10 = 1/bf01;
            
            if nargout>1
                % Traditional tests
                pHat = binofit(k,n,0.05);
                p = 1-binocdf(k,n,p);
            end
            
        end
        
        
        function [bf10,r,p] = corr(arg1,arg2)
            % Calculate the Bayes Factor for Pearson correlation between two
            % variables.
            % INPUT
            % (x,y)  - two vectors of equal length.
            % OR
            % (r,n)  - the correlation and number of samples
            %
            % OUTPUT
            % bf10 = the Bayes Factor for the hypothesis that r is differnt
            %           from zero (two-tailed).
            % r - the correlation
            % p - the tradiational p-value based on Fisher-transformed
            %
            if isscalar(arg1) && isscalar(arg2)
                r= arg1;
                n= arg2;
            else
                x=arg1;y=arg2;
                [r,p] = corr(x,y,'type','pearson');
                n=numel(x);
            end
            
            % Code from Sam Schwarzkopf
            F = @(g,r,n) exp(((n-2)./2).*log(1+g)+(-(n-1)./2).*log(1+(1-r.^2).*g)+(-3./2).*log(g)+-n./(2.*g));
            bf10 = sqrt((n/2)) / gamma(1/2) * integral(@(g) F(g,r,n),0,Inf);
            
            if nargout>1
                % Compute classical stats too
                t = r.*sqrt((n-2)./(1-r.^2));
                p = 2*tcdf(-abs(t),n-2);
            end
        end
    end
    
    %% Hide some of the handle class member functions for ease of use.
    methods (Hidden=true)
        function notify(o)
        end
        function addlistener(o)
        end
        function findobj(o)
        end
        function findprop(o)
        end
        function listener(o)
        end
        
        
        
        
    end
end