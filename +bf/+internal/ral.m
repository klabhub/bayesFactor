classdef ral
    % Class ral - Real As Log
    %
    % Improves numerical summation/products/integration for very small or
    % very large numbers (e.g. products of probabilities) by representing
    % them as sgn(x)*exp(log(abs(x))) and using numerically stable
    % computations for elementary operations (+,-,/,*)
    %
    % Each object represents a single real number.
    %
    % The improved ability to handle low-probability events comes at a cost
    % of speed of computation.
    %
    % BK - Dec 2020
    
    properties (SetAccess={?bf.internal.ral}, GetAccess=public)
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
                x1 = log(abs(x1));
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
        
        % needed for integration with @integral,
        % but that would also require implemnting linear algebra...
        % (mtimes)
        %         function v = isfloat(o)
        %             v = true;
        %         end
        %         function v=superiorfloat(o,x)
        %             v = 'double';
        %         end
        %
        function v = double(o)
            % Convert a ral to a double.
            % Of course this could result in Inf even when the ral is
            % not inf (but has high pwr).
            isz = o.isZero;
            v = nan(size(o));
            v(isz) = 0;
            os= o.sgnArray;
            op = o.pwrArray;
            v(~isz) = os(~isz).*exp(op(~isz));
        end
        
        function [v] = num2str(o,n,asRal)
            % For use in axis labels etc.
            % default output is s*exp(pwr)
            % n determines the precision of the pwr.[3]
            % Set asRal to false to get a standard double representation [true]
            nin =nargin;
            if nin<3
                asRal = true;
                if nin< 2
                n =3;
                end;
            end
            if asRal
                sLabels = {'-','0',' '};
                v = sprintf(['%sexp(%.' num2str(n) 'f)'],sLabels{o.sgn+2},o.pwr);
            else
                v = sprintf(['%.' num2str(n) 'f'],double(o));
            end        
    end
    
    function disp(o)
    %Show the list of numbers in their ral format
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
    
    
    function v =pwrArray(o)
    % Extract an array of pwr values from an array of ral objects
    v = [o.pwr];
    v = reshape(v,size(o));
    end
    function v =sgnArray(o)
    % Extract an array of sgn values from an array of ral objects.
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
        v = copy(o,[],-o.sgnArray);
    end
    
    function v = reciprocal(o)
        v =copy(o,-o.pwrArray);
    end
    
    function v= power(o,e)
        assert(isa(e,'double'),'Exponent must be a double');
        if numel(e) ~=1 && any(size(e) ~=size(o))
            error(['The size of the exponent ' num2str(size(e)) ' does not match the size of the RAL array' num2str(size(o)) ]);
        end
        if numel(e)==1
            e = repmat(e,size(o));
        end
        np = nan(size(o));
        ns = nan(size(o));
        
        op  =o.pwrArray;
        os = o.sgnArray;
        
        zeroExp = e==0;
        np(zeroExp) = 0;
        ns(zeroExp) = 1;
        isEven = floor(e)==round(e) & mod(e,2)==0;
        np(isEven) = op(isEven).*double(e(isEven));
        ns(isEven) = 1;
        rest= isnan(np);
        np(rest) = op(rest).*double(e(rest));
        ns(rest) = os(rest);
        v = copy(o,np,ns);
    end
    
    function v = sqrt(o)
        v = power(o,0.5);
    end
    
    function v = plus(o,x)
        [op,os,xp,xs] = getArrays(o,x);
        [np,ns] = bf.internal.ral.plusRal(op,os,xp,xs);
        v = copy(o,np,ns);
    end
    
    
    function v= minus(o,x)
        [op,os,xp,xs] = getArrays(o,x);
        [np,ns] = bf.internal.ral.minRal(op,os,xp,xs);
        v = copy(o,np,ns);
    end
    
    function v = times(o,x)
        [op,os,xp,xs] = getArrays(o,x);
        np = op+xp;
        ns = os.*xs;
        v = copy(o,np,ns);
    end
    
    function v = mtimes(o,x)
        persistent warned
        if isempty(warned)
            warned =true;
            fprintf(2,'Matrix multiplication for RAL not implemented yet...assuming you want to do .* instead of * (Further warnigns are disabled for this session)');
        end
        v = times(o,x);
    end
    
    function v = rdivide(o,x)
        [op,os,xp,xs] = getArrays(o,x);
        np = op-xp;
        ns = os.*xs;
        v = copy(o,np,ns);
    end
    
    function v = mrdivide(o,x)
        persistent warned
        if isempty(warned)
            warned =true;
            fprintf(2,'Matrix division for RAL not implemented yet...assuming you want to do ./ instead of / (Further warnigns are disabled for this session)');
        end
        v = rdivide(o,x);
    end
    
    %% Relational operators
    function v =isZero(o)
        v = (isinf(o.pwrArray) & sign(o.pwrArray)==-1) | o.sgnArray==0;
    end
    
    function v =eq(o,x)
        [op,os,xp,xs] = getArrays(o,x);
        v = op==xp & os==xs;
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
        [op,os,xp,xs] = getArrays(o,x);
        v = bf.internal.ral.gtPwrSgn(op,os,xp,xs);
    end
    
    
    function [m,s] = mstd(o,dim)
        % Determine mean and standard deviation
        %
        % https://www.johndcook.com/blog/standard_deviation/
        % INPUT
        % dim - dimension along which to operate ([1])
        % NOTE
        % Unlike the standard mean function, this one does not
        % automatically transpose a row vector. (I.e. dim is 1 by
        % default and for a row vector the mean will simply return the
        % entries in each column).
        %
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
                m = oldM + (o(i,allDims{:}) - oldM)./i;
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

%% Internal utility member functions
methods (Access=protected)
    
    function v = copy(o,np,ns)
        % Copying from an existign o array is faster than creating a
        % new one, so we use this to return results from many elmat
        % operations.
        % INPUT
        %  o = array of ral objects
        % np - new power values for each element in o.
        % ns - new sign values ofr each element in o.
        %
        % OUTPUT
        % v = array of ral objects with the new power and sign values.
        %       If np or ns is empty (or omitted), then v will have the
        %       corresponding values from o.
        v=o;
        if nargin >1 && ~isempty(np)
            np = num2cell(np);
            [v.pwr] = deal(np{:});
        end
        if nargin >2 && ~isempty(ns)
            ns = num2cell(ns);
            [v.sgn] = deal(ns{:});
        end
    end
    
    function [op,os,xp,xs] = getArrays(o,x)
        % Used by arithmetic operations to do singleton expansion and
        % allow operations on two rals (instead of a ral and a double).
        
        % Extract or calculate pwr and sgn from o and x
        if bf.internal.ral.isral(o)
            os = o.sgnArray;
            op = o.pwrArray;
        else
            os = sign(o);
            op = log(abs(o));
        end
        if bf.internal.ral.isral(x)
            xs = x.sgnArray;
            xp = x.pwrArray;
        else
            xs = sign(x);
            xp = log(abs(x));
        end
        
        % Singleton expansion (requires more memory but makes
        % subsequent code a lot easier).
        if numel(xs)==1 && ~numel(os)==1
            xs = repmat(xs,size(os));
            xp = repmat(xp,size(os));
        end
        if numel(os)==1 && ~(numel(xs)==1)
            os = repmat(os,size(xs));
            op = repmat(op,size(xs));
        end
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
        % The calling code should make sure that x>y
        assert(all(x>y,'all'),'x>y for log(exp(x)-exp(y))');
        v = x + log(1-exp(y-x));
    end
    
    
    function  v =logExpXPlusExpY(x,y)
        % v = log(exp(x) +exp(y))
        xGty = x>y;
        v = nan(size(x));
        delta= x-y;
        v(xGty) =x(xGty) + bf.internal.ral.log1PlusExpX(-delta(xGty));
        xLty = x<y;
        v(xLty) =y(xLty) + bf.internal.ral.log1PlusExpX(delta(xLty));
        eq = y==x;
        v(eq) = log(2)+x(eq);
    end
    
    function v = gtPwrSgn(op,os,xp,xs)
        % Determine o>x based on arrays of pwr and sgn values.
        % INPUT
        % op = pwr of o
        % os = sgn  of o
        % xp = pwr of x
        % xs = sgn of x
        % OUTPUT
        % v = logical stating o>x
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
    
    function [np,ns] = plusRal(op,os,xp,xs)
        % compute v = o+x
        % Utlity function that works with arrays of pwr and sign
        % (instead of arrays of ral objects). This speeds up the
        % computations considerably as we avoid repeated object
        % creation. See plus(o,x).
        % INPUT
        % op = pwr of o
        % os = sgn  of o
        % xp = pwr of x
        % xs = sgn of x
        % OUTPUT
        % np = pwr of the resulting rals.
        % ns = sgn of the resulting rals.
        
        np=op;
        ns=os;
        oIsZero = os==0;
        xIsZero = xs==0;
        ns(oIsZero) = xs(oIsZero);
        np(oIsZero) = xp(oIsZero);
        
        bothNeg = os ==-1 & xs ==-1;
        if any(bothNeg)
            % v = -log(exp(x)+exp(y))
            np(bothNeg) = bf.internal.ral.logExpXPlusExpY(op(bothNeg),xp(bothNeg));
            ns(bothNeg) = -1;
        end
        
        xNeg = os==1 & xs==-1;
        if any(xNeg)
            % Make both positive then call min
            % o - x -> o + (-x)
            [np(xNeg),ns(xNeg)] = bf.internal.ral.minRal(op(xNeg),os(xNeg),xp(xNeg), -xs(xNeg));
        end
        
        oNeg = os==-1 & xs ==1;
        if any(oNeg)
            % Make both positive then call min
            % o - x -> (-o) + x
            [np(oNeg),ns(oNeg)] =bf.internal.ral.minRal(xp(oNeg),xs(oNeg),op(oNeg),-os(oNeg));
        end
        
        bothPos = ~(oIsZero|xIsZero|bothNeg|xNeg|oNeg);  %
        if any(bothPos)
            % v = log(exp(x) + exp(y))
            np(bothPos) = bf.internal.ral.logExpXPlusExpY(op(bothPos),xp(bothPos));
            ns(bothPos) = 1;
        end
        
        % Force sgn=0 for exp(-inf)
        ns(isinf(np)&np<0) = 0;
    end
    
    
    function [np,ns] = minRal(op,os,xp,xs)
        % compute v = o-x
        % Utlity function that works with arrays of pwr and sign
        % (instead of arrays of ral objects). This speeds up the
        % computations considerably as we avoid repeated object
        % creation. See minus(o,x).
        % INPUT
        % op = pwr of o
        % os = sgn  of o
        % xp = pwr of x
        % xs = sgn of x
        % OUTPUT
        % np = pwr of the resulting rals.
        % ns = sgn of the resulting rals.
        np = op;
        ns = os;
        oIsZero  = os==0;
        xIsZero = xs==0;
        np(oIsZero) = xp(oIsZero);
        ns(oIsZero) = -xs(oIsZero);
        
        bothNeg = os ==-1 & xs ==-1;
        if any(bothNeg)
            % o-x -> -o + (-x)
            [np(bothNeg),ns(bothNeg)] = bf.internal.ral.plusRal(op(bothNeg),-abs(os(bothNeg)),xp(bothNeg),abs(xs(bothNeg)));
        end
        
        xNeg = os==1 & xs==-1;
        if any(xNeg)
            % o-x -> o + (-x)
            [np(xNeg),ns(xNeg)] = bf.internal.ral.plusRal(op(xNeg),os(xNeg),xp(xNeg),-xs(xNeg));
        end
        
        oNeg = os==-1 & xs ==1;
        if any(oNeg)
            % o-x -> -((-o) + x)
            [np(oNeg),ns(oNeg)] = bf.internal.ral.plusRal(op(oNeg),-os(oNeg),xp(oNeg),xs(oNeg));
            ns(oNeg) = -ns(oNeg);
        end
        
        bothPositive = ~(oIsZero|xIsZero|bothNeg|xNeg|oNeg);
        if any(bothPositive)
            %  oGtx = o>x;
            oGtx = bf.internal.ral.gtPwrSgn(op(bothPositive),os(bothPositive),xp(bothPositive),xs(bothPositive));
            ixBothPos =find(bothPositive);
            bothPosAndoGtx = ixBothPos(oGtx);
            np(bothPosAndoGtx) = bf.internal.ral.logExpXMinusExpY(op(bothPosAndoGtx),xp(bothPosAndoGtx));
            ns(bothPosAndoGtx) = 1;
            
            % oLtx = o<x;
            oLtx = bf.internal.ral.gtPwrSgn(xp(bothPositive),xs(bothPositive),op(bothPositive),os(bothPositive));
            
            bothPosAndoLtx = ixBothPos(oLtx);
            np(bothPosAndoLtx) = bf.internal.ral.logExpXMinusExpY(xp(bothPosAndoLtx),op(bothPosAndoLtx));
            ns(bothPosAndoLtx) = -1;
            
            same = ixBothPos(~(oGtx|oLtx));
            np(same) = -inf;
            ns(same) = 0;
        end
        
        % Force sgn=0 for exp(-inf)
        ns(isinf(np)&np<0) = 0;
        
    end
    
    
    %% Static utility functions
    function v = isral(x)
        v = isa(x,'bf.internal.ral');
    end
    
    
    function test
        % Run some tests to see everything is (still) working
        
        
        tol = 10*eps;
        nrRows = 100;
        nrCols = 100;
        r= randn(nrRows,nrCols);
        a= bf.internal.ral(r);
        
        lasterr(''); %#ok<LERR>
        assert(all(a-a< tol,'all'),'Subtraction test Failed');
        assert(all(a+a - 2*r < tol,'all'),'Addition test Failed');
        assert(all(mean((a.*2 - 2*r))<tol),'Mean test failed');
        assert(all(mean(a,1)-mean(r,1) < tol,'all'),'Mean test failed');
        
        if strcmpi(lasterr,'') %#ok<LERR>
            fprintf('All ral tests passed successfully.\n');
        end
    end
end
end