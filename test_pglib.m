clearvars; clc;

% Numerical tolerance for comparisons
tol = 1e-12;

% Add path
pglib_dir = [cd '\pglib'];
addpath(pglib_dir);

% Get list of valid cases from folder
file_list = ls('pglib');
i = size(file_list, 1);
cases = cell(i, 1);
c = 0;
while i > 0
    filename = regexp(file_list(i,:), 'pglib_opf_\w*.m', 'match');
    if isempty(filename)
        % Invalid file, delete from list
        file_list(i,:) = [];
    else
        c = c + 1;
        cases{c} = filename{1}(1:(end-2));
    end
    i = i - 1;
end
cases = cases(c:-1:1);

% Run test cases
imax = size(cases, 1);
flags = zeros(imax, 1);
times = zeros(imax, 1);
nodes = zeros(imax, 1);
lines = zeros(imax, 1);
for i = 1:imax
    case_i = str2func(cases{i});
    mpc = case_i();
    tic;
    [flags(i), nodes(i), lines(i)] = check_inv(mpc, tol); % main function
    times(i) = toc;
end

% Save data
save('results.mat', 'cases', 'flags', 'times', 'nodes', 'lines');

% About systems violating the passivity asumption:
% Only 4 system have inductive branches with negative resistance (1 each).
% All of them have multiple capacitive branches with negative resistance,
% except for 1 system that only has 1 branch with negative resistance (a
% capacitive branch).

% Remove path
rmpath(pglib_dir);

% Print results
print_results();
