function adj_list = bi_adj_list(lines, N)
%BI_ADJ_LIST Bidirectional adjacency matrix of power system
%   adj_list = bi_adj_list(lines, N)
%   Computes the bidirectional adjacency matrix of a network graph with N
%   nodes and L edges described by the L x 4 matrix lines. Parallel
%   branches must be contiguous rows in lines.
%   BI_ADJ_LIST returns a N x 1 cell array adj_list, where the i-th entry 
%   is a structure with the following properties:
%       -edges:
%       Vector whos entries are the neighbor nodes of i.
%       -y:
%       Cell array where each entry is a vector with the admittances of all
%       parallel edges conneting to a specific neighbor. More specifically,
%       adj_list{i}.y{k} is a vector with the admittances of all parallel
%       edges connecting node i with node adj_list{i}.edges{k}. The 
%       admittance values are reflected to the transformer side connected 
%       to node i.
%       -Ybus:
%       Cell array where each entry is the 2x2 admittance matrix formed by
%       all parallel edges connecting node i to a specific neighbor node. 
%       More specifically, adj_list{i}.Ybus{k} the 2x2 admittance matrix 
%       formed by all parallel edges connecting node i with node 
%       adj_list{i}.edges{k}. The first row/column of the matrix 
%       corresponds to node i, and the second row/column corresponds to
%       node adj_list{i}.edges{k}.
%       -rev_ind:
%       Vector with indices of the reverse edges of node i. This indices 
%       are such that adj_list{i}.edge(k) and
%       adj_list{adj_list{i}.edge(k)}.edge(adj_list{i}.rev_ind(k)) describe
%       the same bidirectional edge.
%   This algorithm runs in O(N+L) time.


%% Initialize structure, O(N)
adj_list = cell(N, 1);
for i = 1:N
    adj_list{i}.edges = [];
    adj_list{i}.y = cell(0, 1);
    adj_list{i}.Ybus = cell(0, 1);
    adj_list{i}.rev_ind = [];
end


%% Compute bidirectional adjacency list, O(N+L)
L = size(lines, 1);
for i = 1:L
    y = lines(i,3);
    a = lines(i,4);
    Ybus = [y -a*y; -conj(a)*y (abs(a).^2)*y]; % 2x2 Y bus
    if i > 1 && lines(i,1) == lines(i-1,1) && lines(i,2) == lines(i-1,2)
        % Parallel branch, update edge
        adj_list{lines(i,1)}.y{end}(end+1,1) = y;
        adj_list{lines(i,1)}.Ybus{end} = ...
            adj_list{lines(i,1)}.Ybus{end} + Ybus;
        % Update reverse edge
        adj_list{lines(i,2)}.y{end}(end+1,1) = y ./ (abs(a).^2);
        adj_list{lines(i,2)}.Ybus{end} = ...
            adj_list{lines(i,2)}.Ybus{end} + Ybus([2 1], [2 1]);
    else
        % New edge to be added
        adj_list{lines(i,1)}.edges(end+1) = lines(i,2);
        adj_list{lines(i,1)}.y{end+1} = y;
        adj_list{lines(i,1)}.Ybus{end+1} = Ybus;
        % Reflect admittance and invert taps ratio for the reverse edge.
        adj_list{lines(i,2)}.edges(end+1) = lines(i,1);
        adj_list{lines(i,2)}.y{end+1} = y ./ (abs(a).^2);
        adj_list{lines(i,2)}.Ybus{end+1} = Ybus([2 1], [2 1]);
        % Store indices of reverse edges
        adj_list{lines(i,1)}.rev_ind(end+1) = ...
            length(adj_list{lines(i,2)}.edges);
        adj_list{lines(i,2)}.rev_ind(end+1) = ...
            length(adj_list{lines(i,1)}.edges);
    end
end

