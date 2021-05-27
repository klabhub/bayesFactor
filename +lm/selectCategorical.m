function T= selectCategorical(T,varargin)
% Select rows in a table based on matching  one or more of its categorical 
% columns, and prune the categoricals to remove categories that no longer
% occur.
% 
% INPUT 
% T - A Table
% parm/value pairs identifying rows of the table.
% OUTPUT
% T - A subset of the rows in the table
% 
% BK - May 2021

keep = true(height(T),1);
for i=1:2:numel(varargin)
keep = keep & ismember(T.(varargin{i}),categorical(varargin{i+1}));
if ~any(keep)
    error(['No match with ''' char(varargin{i+1})  ''' for ''' varargin{i} '''. No rows left in the table']);
end
end
T =T(keep,:);

for i=1:2:numel(varargin)
    newCategories = unique(T.(varargin{i}));
    T.(varargin{i}) = setcats(T.(varargin{i}),string(newCategories));
end

end