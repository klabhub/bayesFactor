% Utility script to publish the gettingStarted.m script as index.html on
% GitPages.
options = struct('format','html','ouputDir','./docs/');
htmlDoc =publish('gettingStarted.m',options);
movefile(htmlDoc,'./docs/index.html')