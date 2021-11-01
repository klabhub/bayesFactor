classdef pWaitBar <handle
    properties 
        queue;
        handle;
        total;
        current;
        message = 'Please wait...';
    end
    methods
        function o=pWaitBar(totalNr,message)         
            o.total = totalNr;
            o.current = 0;
            o.message = message;
            o.queue = parallel.pool.DataQueue;
            afterEach(o.queue,@(x) o.update(x));
            o.handle = waitbar(0,o.message);
            update(o);
        end
        function increment(o,x)
            % x is the item # in the increment call. It is ignored because
            % its number is not meaningful in a parfor
            o.current = o.current +1;   
            send(o.queue,x);            
        end
        function update(o,~)
            waitbar(o.current/o.total,o.handle);
            o.current
        end
        function close(o)
            close(o.handle);
        end
        function delete(o)
            close(o);
        end
    end
end
