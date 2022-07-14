function [flag, N, L0, Lr] = check_inv(mpc, tol)
%CHECK_INV Check invertibility of admittance matrix
%   flag = check_inv(mpc, tol)
%   This function attempts to apply Theorem 1 of the paper to certify the
%   invertibility of the admittance matrix of a power system with N nodes
%   and L0 lines, Lr of which are purely reatice. mpc is a struct with all 
%   the power system information, in MATPOWER format. tol is a tolerance 
%   used for numerical comparations (mostly to determine wheter a given 
%   values is zero or not). CHECK_INV returns an integer flag, which takes 
%   one of the following values:
%   SUCCESS VALUES:
%       1:
%       According to Theorem 1 the admittance matrix is invertible.
%       2:
%       According to Theorem 1 the admittance matrix is singular.
%   CORNER-CASE VALUES:
%       0:
%       There are no shunts, so according to Theorem 1 the rank of the
%       admittance matrix equals the rank of the generalized incidence
%       matrix. The matrix may be singular.
%   FAILURE VALUES:
%       -1:
%       Assumption 1 is violated: after reduction some series branches are 
%       zero or non-finite.
%       -2:
%       Assumption 2 is violated: some taps ratio are zero.
%       -4:
%       Passivity assumption of Theorem 2 is violated: after reduction some
%       admittances were found to have negative conductance.
%       -5:
%       Connectivity assumption of Theorem 1 is violated: the network graph
%       is not connected.
%       -6:
%       Some reactive component of the network does not have any of the 
%       required topologies of Theorem 3, the algorithm cannot proceed.
%       -7:
%       Some reactive component of the network satisfies the topology of
%       condition 1) of Theorem 3, but none of the conditions are satisfied
%       anyway, so the algorithm cannot proceed.
%   A failure value indicates that the invertibility of the admittance
%   matrix cannot be asserted. A sucess value indicates that the function
%   sucessfully certified that the matrix is invertible.
%   This function runs in O(N+L) time.


%% Extract data from MATPOWER struct, O(N+L)
Sb = mpc.baseMVA; % base power
N = size(mpc.bus, 1); % system size (# of nodes)
labels = mpc.bus(:,1); % bus labels

% Line data
lines = mpc.branch(:,[1:5 9 10]);
lines(lines(:,6)<=tol,6) = 1; % change zero taps by 1
L = size(lines, 1);
L0 = L;

% Replace bus labels with 1:N numbering in line data
inv_labels = sparse(labels, ones(N, 1), 1:N);
n_labels = length(inv_labels);
A = sparse(1:L, lines(:,1), ones(L, 1), L, n_labels);
lines(:,1) = A * inv_labels; % O(L)
A = sparse(1:L, lines(:,2), ones(L, 1), L, n_labels);
lines(:,2) = A * inv_labels; % O(L)

% Shunt data
shunts = (mpc.bus(:,5) + 1i.*mpc.bus(:,6)) ./ Sb;
for i = 1:L
    % Add line shunts
    shunts(lines(i,1)) = shunts(lines(i,1)) + 0.5i .* ...
        lines(i,5) ./ (lines(i,6).^2);
    shunts(lines(i,2)) = shunts(lines(i,2)) + 0.5i .* lines(i,5);
end

% Rewrite line data using the following cols now:
% lines=[ from_bus  to_bus  series_admittance  complex_taps_ratio];
lines = [lines(:,1:2) 1./(lines(:,3)+1i.*lines(:,4)) ...
    lines(:,6).*exp((1i.*pi./180).*lines(:,7))];

% Count the number of purely reactive branches
Lr = sum(abs(real(lines(:,3))) <= tol);

% Check Assumption 2
if any(abs(lines(:,4)) <= tol)
    flag = -2;
    return;
end


%% Data processing O(N+L)
% Invert start and end nodes to make the end node always larger
for i = 1:L
    if lines(i,1) > lines(i,2)
        % Reflect impedances and admittances
        lines(i,3) = lines(i,3) ./ abs(lines(i,4).^2);
        % Invert taps ratio
        lines(i,4) = 1./lines(i,4);
        % Invert nodes
        temp = lines(i,2);
        lines(i,2) = lines(i,1);
        lines(i,1) = temp;
    end
end

% Sort taps ratio column with radix sort, O(L)
% taps = round(lines(:,4) ./ (2.*tol));
taps = floor(lines(:,4) ./ (2.*tol)); % remove remainder mod 2*tol
lines(:,4) = 2.*tol .* taps;
bus_perm = radix_sort(imag(taps) - min(imag(taps))); % O(L)
lines = lines(bus_perm,:);
taps = taps(bus_perm);
bus_perm = radix_sort(real(taps) - min(real(taps))); % O(L)
lines = lines(bus_perm,:);

% Use counting sorting to sort the end node column
bus_perm = counting_sort(lines(:,2), N); % O(N+L)
lines = lines(bus_perm,:);

% Now we use counting sort (which is stable) for the start node column
bus_perm = counting_sort(lines(:,1), N); % O(N+L)
lines = lines(bus_perm,:);

% Combine parallel branches, if possible
i = 1;
while i < L
    if lines(i,1) == lines(i+1,1) && lines(i,2) == lines(i+1,2)
        % Parallel branches found
        if abs(lines(i,4) - lines(i+1,4)) <= tol
            % Combine branches and delete one of them
            lines(i,3) = lines(i,3) + lines(i+1,3);
            lines(i+1,:) = [];
            L = L - 1;
        else
            % Parallel branches cannot be reduced
            i = i + 1;
        end
    else
        i = i + 1;
    end
end

% Check Assumption 1
if any(abs(lines(:,3)) <= tol) || any(~isfinite(lines(:,3)))
    flag = -1;
    return;
end

% Check passivity assumption
if any(real(lines(:,3)) < -tol) || any(real(shunts) < -tol)
    flag = -4;
    return;
end

% Compute bidirectional adjacency list, O(N+L)
adj_list = bi_adj_list(lines, N);


%% Find connected components using BFS, O(N+L)
[comp_nodes, trees] = get_comps(adj_list, 1); % Run BFS, O(N+L)

% Check Connectivity
if length(comp_nodes) ~= 1
    flag = -5;
    return;
end


%% Find reactive connected components using BFS, O(N+L)
% Reduce remove lines and shunts with resistive component
ind_reac = abs(real(lines(:,3))) <= tol;
lines_reac = lines(ind_reac,:);
shunts_reac = shunts;
ind_res = abs(real(shunts_reac)) > tol;
shunts_reac(ind_res) = 0;

% Compute bidirectional adjacency list, O(N+L)
adj_reac = bi_adj_list(lines_reac, N);

% Run BFS, O(N+L)
[comp_reac, trees_reac] = get_comps(adj_reac, 1);


%% Check Theorem 3 on each reactive component, O(N+L)
% O(N_k+L_k) per component, O(N+L) in total
c = length(comp_reac);
for k = 1:c
    % Check Condition 3), O(N_k+L_k)
    shunts_k = shunts_reac(comp_reac{k});
    y_k = shunts_k;
    for i = (comp_reac{k}(:).')
        % Branch admittances are included twice (the original and the
        % reflection over the other node), but it does not affect the
        % evaluation of the condition.
        y_k = [y_k; cell2mat(adj_reac{i}.y(:))]; %#ok<AGROW>
    end
    if all(imag(y_k) >= -tol) || all(imag(y_k) <= tol)
        % There are only capacitors or inductors, Condition 3) applies
        continue;
    end
    
    % Check condition 2), O(N_k)
    trees_k = trees_reac{k};
    if all(abs(shunts_k) <= tol) && ~isempty(trees_k)
        % We have a tree with no shunts, Condition 2) applies
        continue;
    end
    
    % Check condition 1), O(N_k+L_k)
    if ~isempty(trees_k)
        % We have a tree
        [cond1_true, fail_node] = ...
            condition1(trees_k, adj_reac, shunts_reac, tol);
        if cond1_true
            continue;
        elseif fail_node ~= trees_k{1}.node
            % Condition 1) is violated for at a node different than root.
            % We can rebuild the tree using fail_node as the new root, if  
            % the condition is still violated then no choice of root  
            % will satisfy the condition.
            % We remark that if the condition is violated at the root node,
            % then no choice of root will satisfy the condition either.
            [~, new_trees] = get_comps(adj_reac, fail_node);
            cond1_true = ...
                condition1(new_trees{1}, adj_reac, shunts_reac, tol);
            if cond1_true
                continue;
            else
                % Some equivalent admittance violates condition 1)
                flag = -7;
                return;
            end
        end
    end
    
    % No condition is satisfied, theorem 3 cannot be applied to this
    % component.
    flag = -6;
    return;
end


%% Apply Theorem 1, O(N)
if any(abs(shunts) > tol)
    % If there are shunts, then the admittance matrix is invertible
    flag = 1;
else
    % There are no shunts
    if ~isempty(trees{1})
        % The network is radial, so the generalized incidence matrix has 
        % N-1 rows and thus the admittance matrix has rank N-1. In
        % particular, the admittance matrix is singular.
        flag = 2;
    else
        % The network is meshed, in this case the code does not decide
        % whether the admittance matrix is invertible or not.
        flag = 0;
    end
end

% % Compute admittance matrix
% ind_i = [lines(:,1); lines(:,2); lines(:,1); lines(:,2); (1:N).'];
% ind_j = [lines(:,2); lines(:,1); lines(:,1); lines(:,2); (1:N).'];
% entries = [-lines(:,3)./conj(lines(:,4)); -lines(:,3)./lines(:,4); 
%     lines(:,3)./abs(lines(:,4).^2); lines(:,3); shunts];
% Y = sparse(ind_i, ind_j, entries, N, N);

