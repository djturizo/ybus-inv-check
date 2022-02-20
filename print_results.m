function print_results(varargin)
%PRINT_RESULTS Print results of invertibility checks in console
%   print_results()
%   print_results(ignore_flags)
%   This function loads the results stored in 'results.mat' and prints them
%   in console. Optional argument ignore_flags is a vector whose entries 
%   are flag values of the function CHECK_INV that will NOT be printed. 
%   This means that test cases with one of the flag values of ignore_flags 
%   will not be printed. When no input is provided no flags are ignored,
%   and all test cases are printed.


%% Load results and filter out test cases to be ignored
% Load results
data = load('results.mat');
flags = data.flags;

% Get indices of test cases to be printed
if nargin > 0
    ignore_flags = varargin{1};
    ind_print = all(bsxfun(@ne, flags, ignore_flags), 2);
else
    ind_print = true(size(flags));
end

% Filter test cases
cases = data.cases(ind_print);
nodes = data.nodes(ind_print);
lines = data.lines(ind_print);
flags = data.flags(ind_print);
times = data.times(ind_print);

% Sort test case by number of nodes
[~, ind] = sort(nodes, 'ascend');
cases = cases(ind);
nodes = nodes(ind);
lines = lines(ind);
flags = flags(ind);
times = times(ind);


%% Print results as a table
fprintf(['_________________________________________________' ...
    '_______________________________\n']);
fprintf(['|                                |       |       ' ...
    '|Thm.   |Is YN |     |        |\n']);
fprintf(['|Test case                       |N      |L      ' ...
    '|holds? |inv.? |Flag |Time[s] |\n']);
fprintf(['+--------------------------------+-------+-------' ...
    '+-------+------+-----+--------+\n']);
sep = '|'; % column separator
% sep = '&'; % column separator
kmax = length(cases);
for k = 1:kmax
    fprintf(sep);
    fprintf(['%-31s ' sep '%-6d ' sep '%-6d '], ...
        cases{k}(1:min(end,31)), nodes(k), lines(k));
    if flags(k) > 0
        fprintf([sep 'Yes    ']);
        if flags(k) == 0
            fprintf([sep '?     ']);
        else
            fprintf([sep 'Yes   ']);
        end
    else
        fprintf([sep 'No     ' sep '-     ']);
    end
    fprintf([sep '%-4d ' sep '%-7.2f ' sep], flags(k), times(k));
    fprintf('\n');
%     fprintf('\\\\\n\\hline\n');
end
fprintf(['+--------------------------------+-------+-------' ...
    '+-------+------+-----+--------+\n']);
fprintf('\n');

