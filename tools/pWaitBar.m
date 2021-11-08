classdef pWaitbar < parallel.pool.DataQueue
    %NOT FUNCTIONAL - 
    % For some reason this does not work... 
    % Dont know why, but run the example code below to get the error
    % messages.
    % 
    % A waitbar that works with parralel processing using parfor.
    % USAGE
    % Create the object before the parfor loop using the same format as
    % waitbar 
    % h =  pWaitBar(100,'Please wait...');
    % Then let the workers call 
    % increment(h)
    % after each loop through the parfor block.
    % 
    % Example:
    %     h=pWaitBar(100,'Please wait');
    %     parfor i=1:100
    %        m=magic(5); % Do something
    %        increment(h,1); % Another one done.
    %     end
    %     close(h)
    %     % BK - Nov 21
    properties 
        handle;
        total;
        current;
        message = 'Please wait...';
    end
    methods
        function o=pWaitbar(totalNr,message)         
            % Create a waitbar for use with parfor
            % INPUT
            % totalNr 
            % message 
            % OUTPUT 
            % o = object.
            error('pWaitbar is not functional yet...')
            if nargin<2
                message = 'Please wait...';
                if nargin <1
                    totalNr = 100;
                end
            end
               
            o.total = totalNr;
            o.current = 0;
            o.message = message;
            afterEach(o,@(x) o.update(x));
            o.handle = waitbar(0,o.message);
            update(o);
        end
        function increment(o,x)
            % Increment the interal counter (1 more job done) and update
            % the waitbar.
            % 
            % x is the optional item # in the increment call. It is ignored because
            % its number is not meaningful in a parfor
            if nargin <2 
                x = 1;
            end
            o.current = o.current +1;   
            o.send(1);            
        end
        function update(o,~)
            waitbar(o.current/o.total,o.handle);
            o.current
        end
        function close(o)
            if ishandle(o.handle)
                close(o.handle);
            end
        end
        function delete(o)
            close(o);
        end
    end
end
