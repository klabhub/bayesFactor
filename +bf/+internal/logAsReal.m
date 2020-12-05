classdef logAsReal
    properties
        sgn; % Sign of the number.
        pwr; %
    end
    
    methods
        function v= numel(o)
            v =numel(o.pwr);
        end
        function v =size(o,dim)
            if nargin >1
                v= size(o.pwr,dim);
            else
                v = size(o.pwr);
            end
        end
        
        function v = isnan(v)
            v = isnan(o.pwr);
        end
        function disp(o)
            for i=1:numel(o)
                fprintf('%d*exp(%3.3g)\n',o.sgn(i),o.pwr(i));
            end
        end
        
        function o = logAsReal(p,s)
            if nargin <2
                % No signs provided, interpret as a real
                o = bf.internal.logAsReal(log(abs(p)),sign(p));
            else
                % Both pwr and sign provided
                if ~all(ismember(s(~isnan(s)),[-1 0 1]))
                    error('Sign must be -1 0 or 1 not %g',s(~isnan(s)));
                end
                expMinusInf = isinf(p) & sign(p)==-1;% exp(-inf)
                s(expMinusInf) = 0;
                p(s==0)= -Inf;
                o.sgn = s;
                o.pwr= p;
            end
        end
        
        function v= log(o)
            switch (o.sgn)
                case 1
                    v = o.pwr;
                case 0
                    v = -inf;
                case -1
                    v = nan;
                otherwise
                    error('Sign must be -1 0 or 1 not %g',o.sgn);
            end
            
        end
        
        function v = abs(o)
            v= bf.internal.logAsReal(o.pwr,1);
        end
        
        
        function v = uminus(o)
            v = o.negative;
        end
        function v =negative(o)
            v= bf.internal.logAsReal(o.pwr,-o.sgn);
        end
        
        function v = reciprocal(o)
            v = bf.internal.logAsReal(-o.pwr,o.sgn);
        end
        
        
        function v= ref(o,ix)
            if numel(o)==1 && islogical(ix)
                ix=true;
            end
            v = subsref(o,substruct('()',ix));
        end
        function o= asgn(o,ix,v)
            if numel(o)==1 && islogical(ix)
                ix=true;
            end 
            o= subsasn(o,substruct('()',ix,v));
        end
        function v = subsref(o,subs)
            switch (subs.type)
                case '()'
                    p = o.pwr(subs.subs);
                    s = o.sgn(subs.subs);
                    v= bf.internal.logAsReal(p,s);
                otherwise
                    error('');
            end
        end
        
        function v = subsasgn(o,subs,v)
            switch (subs.type)
                case '()'
                    o.pwr(subs.subs) =v.pwr;
                    o.sgn(subs.subs) =v.sgn;                    
                otherwise
                    error('');
            end
        end
        
        
        
        function v= power(o,e)
            assert(numel(e)==1,'Exponentiation not vectorized');
            assert(isa(e,'double'),'Exponentiation only with doubles')
            p = nan(size(o));
            s = nan(size(o));
            zeroExp = e==0;
            p(zeroExp) = 0;
            s(zeroExp) = 1;
            isEven = floor(e)==round(e) & mod(e,2)==0;
            p(isEven) = o.pwr(isEven).*e(isEven);
            s(isEven) = 1;
            rest= isnan(p);
            p(rest) = o.pwr(rest).*e(rest);
            s(rest) = o.sgn(rest);
            v =bf.internal.logAsReal(p,s);
        end
        function v =isZero(o)
            v = (isinf(o.pwr) & sign(o.pwr)==-1) | o.sgn==0;
        end
        
        function v = double(o)
            v = o.sgn.*exp(o.pwr);
        end
        
        function v= select(o,x)
            v = bf.internal.logAsReal(o.pwr(x),o.sgn(x));
        end
        function v = plus(o,x)
            if isa(x,'bf.internal.logAsReal')
                
                v = bf.internal.logAsReal(nan(size(o)));
                v(o.isZero) = x(o.isZero);
                v(x.isZero) = o(x.isZero);
                
                bothNeg = o.sgn ==-1 & x.sgn ==-1;
                if any(bothNeg)
                    if numel(x)==1
                        tmp = ref(o,bothNeg).negative + x.negative;
                    else
                        tmp = ref(o,bothNeg).negative + ref(x,bothNeg).negative;   
                    end
                    v = asgn(v,bothNeg,tmp.negative);
                end
                
                xNeg = o.sgn==1 & x.sgn==-1;
                if any(xNeg)
                    if numel(x)==1
                            v= asgn(v,xNeg, ref(o,xNeg) - x.negative);
                    else
                        v= asgn(v,xNeg, ref(o,xNeg) - ref(x,xNeg).negative);
                    end
                end
                
                oNeg = o.sgn==-1 & x.sgn ==1;                
                if any(oNeg)
                    if numel(x)==1
                        v=asgn(v,oNeg,x - ref(o,oNeg).negative);
                    else
                        v=asgn(v,oNeg, ref(x,oNeg) - ref(o,Neg).negative);
                    end                        
                end
                rest = isnan(v);
                if any(rest)
                    if numel(x)==1                        
                        p = bf.internal.logAsReal.logExpXPlusExpY(o(rest).pwr,x.pwr);
                    else
                        p = bf.internal.logAsReal.logExpXPlusExpY(o(rest).pwr,x(rest).pwr);
                    end
                    v =asgn(v,rest, bf.internal.logAsRal(p,ones(size(p))));
                end
            else
                x = bf.internal.logAsReal(log(abs(x)),sign(x));
                v = o+x;
            end
        end
        
        
        
        function v= minus(o,x)
            if isa(x,'bf.internal.logAsReal')
                
                v = bf.internal.logAsReal(nan(size(o)));
                if any(o.isZero)
                v(o.isZero) = -x(o.isZero);
                end
                if any(x.isZero)
                v(x.isZero) = o(x.isZero);
                end                
                bothNeg = o.sgn ==-1 & x.sgn ==-1;
                if any(bothNeg) %  v = -abs(o) + abs(x);
                    if numel(x)==1
                        tmp = -ref(o,bothNeg).abs + x.abs;
                    else
                        tmp = -ref(o,bothNeg).abs + ref(x,bothNeg).abs;   
                    end
                    v = asgn(v,bothNeg,tmp);
                end
                
                xNeg = o.sgn==1 & x.sgn==-1;                
                if any(xNeg) %v = o+x.negative;
                    if numel(x)==1
                            v= asgn(v,xNeg, ref(o,xNeg) + x.negative);
                    else
                        v= asgn(v,xNeg, ref(o,xNeg) + ref(x,xNeg).negative);
                    end
                end
                
                oNeg = o.sgn==-1 & x.sgn ==1;                
                if any(oNeg) % v = negative(xs + o.negative);
                    if numel(x)==1
                        tmp  = ref(o,oNeg).negative + x;
                    else
                        tmp  = ref(o,oNeg).negative + ref(x,oNeg);
                    end    
                    v=asgn(v,oNeg,tmp.negative);
                end
                
 
                rest = isnan(v);
                if numel(x)==1
                    oGtx = ref(o,rest)>x;
                    oLtx = ref(o,rest)<x;
                else
                    oGtx = ref(o,rest)>ref(x,rest);
                    oLtx = ref(o,rest)>ref(x,rest);
                end
                
                if any(rest & oGtx)   
                    v = asgn(v,rest&oGtx,bf.internal.logAsReal(bf.internal.logAsReal.logExpXMinusExpY(ref(o,rest&oGtx).pwr,ref(x,rest&oGtx).pwr),1));
                end
                if any(rest & oLtx)
                    v = asgn(v,rest&oLtx,bf.internal.logAsReal(bf.internal.logAsReal.logExpXMinusExpY(ref(x,rest&oLtx).pwr,ref(o,rest&oLtx).pwr),-1));
                end
                    v = asgn(v,rest & ~(oGtx|oLtx), bf.internal.logAsReal(zeros(1,sum(rest & ~(oGtx|oLtx)))));
            else
                v = o-bf.internal.logAsReal(x);
            end
        end
        
        function v = lt(o,x)
            % Less than
            if isa(x,'bf.internal.logAsReal')
                v = gt(x,o);
            else
                v = gt(bf.internal.logAsReal(abs(x),sign(x)),o);
            end
        end
        
        function v = gt(o,x)
            if isa(x,'bf.internal.logAsReal')
                if o.sgn ~= x.sgn
                    v = o.sgn > x.sgn;
                elseif o.sgn>0
                    v = o.pwr>x.pwr;
                else
                    v = o.pwr < x.pwr;
                end
            else
                v = gt(o,bf.internal.logAsReal(abs(x),sign(x)));
            end
        end
        
        
        function [m,s] = mstd(o)
            % Determine mean and standard deviation
            %https://www.johndcook.com/blog/standard_deviation/
            m = bf.internal.logAsReal(o(1),1);
            s = bf.internal.logAsReal(0,0);
            for  i = 2:numel(o)
                oldM = m;
                if isnan(o(i))
                    break;
                end
                thisX = bf.internal.logAsReal( o(i), 1);
                m = oldM + ( thisX - oldM) / i;
                s = s + ( thisX - oldM ) * ( thisX - m )/(i-1);
            end
        end
    end
    
    methods (Static)
        
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
            v = x + log(1-exp(y-x));
        end
        
        
        function  v =logExpXPlusExpY(x,y)
            % v = log(exp(x) +exp(y))
            v=x+bf.internal.logAsReal.log1PlusExpX(y-x);
        end
    end
end