classdef ral
    % Class ral - Real As Log
    % Improves numerical summation/products/integration for very small or
    % very large numbers by representing them as sgn(x)*exp(log(abs(x)))
    % Each object represents a single real number, many standard operations
    % have been overloaded to allow addition/substraction, etc.
    %
    % BK - Dec 2020
    
    properties
        sgn; % Sign of the number. (-1 0 1 or NaN)
        pwr; % Power of the number
    end
    
    methods (Access=public)
        
        function o = ral(x1,x2)
            % Constructor for a RAL (Real As Log)
            % x1 = an array of numbers to represent as RALs
            % or
            % x1 : the log of the absolute value of the numbers to
            % represent
            % x2 : the sign of the numbers to represent.
            % OUTPUT
            % An array of RAL objects, matching the size of x1.
            nin=nargin;             
            assert(nin>=0 && nin<3, 'RAL needs 0,1, or 2 inputs , not %d',nin);
            
            if nin==0
                o.sgn =[];
                o.pwr =[];
                return;
            elseif nin==1
                % Single vector of inputs : interpret as a real
                x2 = sign(x1);
                x1 =  log(abs(x1));
            end
            
                % Both pwr and sign provided
                if isempty(x1)
                    o = bf.internal.ral.empty;
                else
                    expMinusInf = isinf(x1) & sign(x1)==-1;% exp(-inf)
                    x2(expMinusInf) = 0;
                    x1(x2==0)= -Inf;
                    ix = 1:numel(x1);
                    s= num2cell(x2);
                    p= num2cell(x1);
                    [o(ix).sgn]  = deal(s{:});                    
                    [o(ix).pwr]  = deal(p{:});
                    o=reshape(o,size(x1));
                end
           
        end
        
        function v = isnan(o)
            v = isnan(o.pwrArray);
        end
        
        function v = double(o)
            isz = o.isZero;
            v = nan(size(o));
            v(isz) = 0;
            os= o.sgnArray;
            op = o.pwrArray;
            v(~isz) = os(~isz).*exp(op(~isz));
        end
        
        function disp(o)
            s= o.sgnArray;
            p = o.pwrArray;
            sLabels = {'-','0',' '};
            %  MAXSHOW= 10;
            disp(['RAL Object with size [',num2str(size(o)) ']'])
            for k=1:size(s,3)
                for i=1:size(s,1)
                    for j=1:size(s,2)
                        if ~isnan(s(i,j,k))
                            fprintf('%sexp(%3.3f)\t',sLabels{s(i,j,k)+2},p(i,j,k));
                        else
                            fprintf('NaN\t')
                        end
                    end
                    fprintf('\n')
                end
            end
        end
    end
    
    methods
        function v =pwrArray(o)
            v = [o.pwr];
            v = reshape(v,size(o));
        end
        function v =sgnArray(o)
            v = [o.sgn];
            v = reshape(v,size(o));
        end
        
    end
    
    
    methods
        %% Elementray Maths
        function [v,asRal]= log(o)
            s= o.sgnArray;
            v = o.pwrArray;
            v(s<0) = NaN;
            v(s==0) = -inf;
            if nargout==2
                asRal = bf.internal.ral(v);
            end
        end
        
        function v = abs(o)
            v= o;
            [v.sgn] = deal(1);
        end
        
        function v = uminus(o)
            v = o.negative;
        end
        function v =negative(o)
            v =o;
            ns = num2cell(-o.sgnArray);
            [v.sgn] = deal(ns{:});            
        end
        
        function v = reciprocal(o)
            v =o;
            np = num2cell(-o.pwrArray);
            [v.pwr] = deal(np{:});              
        end
        
        function v= power(o,e)
            assert(isa(e,'double'),'Exponent must be a double');
            if numel(e) ~=1 && any(size(e) ~=size(o))
                error(['The size of the exponent ' num2str(size(e)) ' does not match the size of the RAL array' num2str(size(o)) ]);
            end
            if numel(e)==1
                e = repmat(e,size(o));
            end
            p = nan(size(o));
            s = nan(size(o));
            
            op  =o.pwrArray;
            os = o.sgnArray;
            
            zeroExp = e==0;
            p(zeroExp) = 0;
            s(zeroExp) = 1;
            isEven = floor(e)==round(e) & mod(e,2)==0;
            p(isEven) = op(isEven).*double(e(isEven));
            s(isEven) = 1;
            rest= isnan(p);
            p(rest) = op(rest).*double(e(rest));
            s(rest) = os(rest);
            np = num2cell(p);
            ns = num2cell(s);
            v = o;
            [v.pwr] = deal(np{:});
            [v.sgn] =deal(ns{:});
        end
        
        
        function v = plus(o,x)
            [o,x] = bf.internal.ral.matchup(o,x);
            os = o.sgnArray;
            xs = x.sgnArray;
            
            v = o;
            v(o.isZero) = x(o.isZero);
            
            bothNeg = os ==-1 & xs ==-1;
            if any(bothNeg)
                v(bothNeg) = -(o(bothNeg).negative + x(bothNeg).negative);
            end
            
            xNeg = os==1 & xs==-1;
            if any(xNeg)
                v(xNeg) = o(xNeg) - x(xNeg).negative;
            end
            
            oNeg = os==-1 & xs ==1;
            if any(oNeg)
                v(oNeg)=x(oNeg) - o(oNeg).negative;
            end
            
            rest = ~(o.isZero|bothNeg|xNeg|oNeg);
         
            if any(rest)
                p = num2cell(bf.internal.ral.logExpXPlusExpY(o(rest).pwrArray,x(rest).pwrArray));
                [v(rest).sgn] = deal(1);
                [v(rest).pwr] = deal(p{:});
            end
        end
        
        
        
        function v= minus(o,x)
            [o,x] = bf.internal.ral.matchup(o,x);
            
            os = o.sgnArray;
            xs = x.sgnArray;
            
            v = o;
            v(o.isZero) = x(o.isZero).negative;
            
            bothNeg = os ==-1 & xs ==-1;
            if any(bothNeg)
                v(bothNeg) = o(bothNeg).abs.negative + x(bothNeg).abs;
            end
            
            xNeg = os==1 & xs==-1;
            if any(xNeg)
                v(xNeg)  = o(xNeg) + x(xNeg).negative;
            end
            
            oNeg = os==-1 & xs ==1;
            if any(oNeg)
                v(oNeg)= -(o(oNeg).negative + x(oNeg));
            end
            
            rest = ~(o.isZero|x.isZero|bothNeg|xNeg|oNeg);
            if any(rest)
                oGtx = o>x;
                restAndoGtx = rest & oGtx;
                p = num2cell(bf.internal.ral.logExpXMinusExpY(o(restAndoGtx).pwrArray,x(restAndoGtx).pwrArray));
                
                [v(restAndoGtx).pwr] = deal(p{:});
                [v(restAndoGtx).sgn] = deal(1);
                oLtx = o<x;
                restAndoLtx = rest & oLtx;
                p = num2cell(bf.internal.ral.logExpXMinusExpY(x(restAndoLtx).pwrArray,o(restAndoLtx).pwrArray));
               [ v(restAndoLtx).pwr] = deal(p{:});
               [ v(restAndoLtx).sgn] = deal(-1);
               
                restSame = rest & ~(oGtx|oLtx);
                
                [v(restSame).pwr]  = deal(-inf);
                [v(restSame).sgn]  = deal(0);
                
            end
            
        end
        
        function v = times(o,x)
            [o,x] = bf.internal.ral.matchup(o,x);
            p = o.pwrArray+x.pwrArray;
            s = o.sgnArray.*x.sgnArray;
            v = bf.internal.ral(p,s);
        end
        
        
        function v = mtimes(o,x)
            persistent warned
            if isempty(warned)
                warned =true;
                fprintf(2,'Matrix multiplication for RAL not implemented yet...assuming you want to do .* instead of *');
            end
            v = times(o,x);
        end
        
        function v = rdivide(o,x)
            [o,x] = bf.internal.ral.matchup(o,x);
            p = o.pwrArray-x.pwrArray;
            s = o.sgnArray.*x.sgnArray;
            v = bf.internal.ral(p,s);
        end
        
        function v = mrdivide(o,x)
            persistent warned
            if isempty(warned)
                warned =true;
                fprintf(2,'Matrix division for RAL not implemented yet...assuming you want to do .* instead of *');
            end
            v = rdivide(o,x);
        end
        
        %% Relational operators
        function v =isZero(o)
            v = (isinf(o.pwrArray) & sign(o.pwrArray)==-1) | o.sgnArray==0;
        end
        
        function v =eq(o,x)
            [o,x] = bf.internal.ral.matchup(o,x);
            v = o.pwrArray==x.pwrArray & o.sgnArray==x.sgnArray;
        end
        
        
        function v = lt(o,x)
            % o<x
            % Compute x>o
            v = gt(x,o);
        end
        
        function v = le(o,x)
            % o<=x
            % Compute ~(o>x)
            v = ~gt(o,x);
        end
        
        function v = ge(o,x)
            % o>=x
            % Compute ~(o<x)
            v = ~lt(o,x);
        end
        
        function v = gt(o,x)
            % o>x
            [o,x] = bf.internal.ral.matchup(o,x);
            v=nan(size(o));
            os = o.sgnArray;
            xs = x.sgnArray;
            op = o.pwrArray;
            xp = x.pwrArray;
            diffSign = os ~=xs;
            v(diffSign) = os(diffSign)>xs(diffSign);
            sameSignPos = ~diffSign & os>0;
            v(sameSignPos) = op(sameSignPos)>xp(sameSignPos);
            sameSignNeg =~diffSign & os<0;
            v(sameSignNeg) = op(sameSignNeg) < xp(sameSignNeg);
            sameSignZero = ~diffSign & os==0;
            v(sameSignZero)  = xs(sameSignZero)<0;            
            v=logical(v);
        end
        
        
        function [m,s] = mstd(o,dim)
            % Determine mean and standard deviation
            %https://www.johndcook.com/blog/standard_deviation/
            nout = nargout;
            if nargin<2
                dim=1;
            end
            nrDims = ndims(o);            
            %if nrDims >= 2 % if 1 or two dimension, don't permute, the output will loose the corect dimensions.
            dims =1:nrDims;
            if dim<=nrDims
                dims(dim) = [];
                dims = [dim dims];
                o= permute(o,dims); % Now the dim-dimension is the first.             
            end
            %end
            allDims = repmat({':'},[1 nrDims-1]);
            nrElms = size(o);
            
            if nrElms(1) ==0
                m=nan;s=nan;
                return
            else
                % Starting point 
                m = o(1,allDims{:});
                if nout>1
                    s = repmat(bf.internal.ral(0,0),[1 nrElms(2:end)]);
                end
                for  i = 2:nrElms(1) 
                    % Loop over the first dimension
                    oldM =m;
                    m = oldM + (o(i,allDims{:}) - oldM) ./ i;
                    if nout>1
                        s = s + (o(i,allDims{:}) - oldM ) * ( o(i,allDims{:}) - m );
                    end                    
                end
                if nout>1
                    s = s/(nrElms(1)-1);
                end
            end
        end
        
        function v = mean(o,dim)
            if nargin <2
                dim =1;
            end
            v = mstd(o,dim);
        end
        function v = std(o,dim)
            if nargin<2
                dim =1;
            end
            [~,v] = mstd(o,dim);
        end
    end
    
    methods (Static)
        %% Static functions to improve numerical stability
        function v= log1PlusExpX(x)
            % v = log(1+exp(x))
            LOGEPS = log(eps('double'));
            LOGQRT = log(0.25);
            if (x > -LOGEPS)
                v =x;
            elseif x>LOGQRT
                v =log(1+exp(x));
            else
                v= log1p(exp(x));
            end
        end
        
        
        function  v = logExpXMinusExpY(x,y)
            % v = log(exp(x) - exp(y))
            %   for x>y only
            assert(all(x>y),'x>y for this formula');
            v = x + log(1-exp(y-x));
        end
        
        
        function  v =logExpXPlusExpY(x,y)
            % v = log(exp(x) +exp(y))
            % 
            xGty = x>y;
            v = nan(size(x));
            delta= x-y;
            v(xGty) =x(xGty) + bf.internal.ral.log1PlusExpX(-delta(xGty));
            xLty = x<y;
            v(xLty) =y(xLty) + bf.internal.ral.log1PlusExpX(delta(xLty));
            eq = y==x;
            v(eq) = log(2)+x(eq);            
        end
        
        %% Static utility functions
        function v = isral(x)
            v = isa(x,'bf.internal.ral');
        end
        
        function [o,x] = matchup(o,x)
            % Make sure both o and x are RAL and
            % do singleton expansion to match the other.
            if ~bf.internal.ral.isral(o)
                o = bf.internal.ral(o);
            end
            if ~bf.internal.ral.isral(x)
                x = bf.internal.ral(x);
            end
            if numel(x)==1
                x = repmat(x,size(o));
            end
            if numel(o)==1
                o = repmat(o,size(x));
            end
        end
        
        function test
            tol = 10*eps;
            nrRows = 100;
            nrCols = 100;
            r= randn(nrRows,nrCols);
            a= bf.internal.ral(r);
            
            assert(all(a-a< tol,'all'),'Subtraction test Failed');
            assert(all(a+a - 2*r < tol,'all'),'Addition test Failed');
            assert(all(mean((a*2 - 2*r))<tol),'Mean test failed');
            assert(all(mean(a)-mean(r) < tol,'all'),'Mean test failed');
            
        end
    end
end