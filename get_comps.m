function [comp_nodes, trees] = get_comps(adj_list, source)
%GET_COMPS Get connected components of a graph
%   comp_nodes = get_comps(adj_list, source)
%   Get the connected components of a graph with N nodes and L edges. 
%   adj_list is the adjacency list of the graph, a N x 1 cell array. The 
%   i-th entry of adj_list is a vector whose entries indicate the nodes 
%   connected to node i. source is an integer indicating the starting node
%   of the BFS algorithm. If there are multiple connected components, new
%   sources will be tried in increasing order.
%   GET_COMPS return a cell array comp_nodes, where each entry containts a 
%   vector with the nodes of a single connected component. Isolated nodes 
%   are not considered as connected components, and are not included in 
%   comp_nodes. 
%   GET_COMPS also returns the cell array trees of the same size as 
%   comp_nodes. The i-th entry of trees is empty is i-th component of the 
%   graph is NOT a tree, otherwise the entry is another cell array with the
%   info of the BFS tree. The k-th entry of the second cell array contains 
%   the info of all nodes which are at distance k-1 from the root of the 
%   BFS tree. More specifically, the k-th entry of the i-th entry of trees 
%   is a struct array with the following fields:
%       -node:
%       Index of the tree node whose info is stored in this entry.
%       -children:
%       Vector with the indices of the children of this node. That is, if
%       trees{i}{k}(h).children(l) = j then trees{i}{k+1}(j).node is a
%       child of trees{i}{k}(h).node in the BFS tree. Notice that from the
%       way the BFS tree is built then all the children of 
%       trees{i}{k}(h).node are in trees{i}{k+1}.
%       -par_ind:
%       Index of the parent of node. More specifically,
%       trees{i}{k-1}(par_ind).node is the parent of node in the BFS tree.
%       -par_edge:
%       Index of the edge connecting node to its parent. More specifically,
%       let p = trees{i}{k-1}(par_ind).node be the parent of node, then
%       adj_list{p}.edges(par_edge) = node.
%   This algorithm runs in O(N+L) time.


%% Initialize variables, O(N)
N = length(adj_list);
visited = false(N, 1);
source_list = 1:N; % List of nodes to use as start of the BFS
source_list(source) = [];
source_list = [source source_list]; % Make source the first on the list
source_ind = 1;
comp_nodes = cell(0, 1);
trees = cell(0, 1);
c = 0; % Number of connected components

% Create list to store unvisited edges (0 if unvisited, 1 if else).
adj_visited = cell(N, 1);
for i = 1:N
    adj_visited{i} = zeros(size(adj_list{i}.edges));
end


%% Run BFS on all components, O(N+L)
while source_ind <= N
    % Create BFS with root on source node
    visited(source) = true;
    c = c + 1;
    comp_nodes{c} = source;
    trees{c} = cell(1, 1);
    trees{c}{1}(1).node = source;
    trees{c}{1}(1).children = [];
    trees{c}{1}(1).par_ind = [];
    trees{c}{1}(1).par_edge = [];
    
    % Create BFS queue with the following info:
    %       [node   depth  par_ind]
    queue = [source 0      1      ];
    % Start BFS
    while ~isempty(queue)
        node = queue(1,1);
        depth = queue(1,2);
        par_ind = queue(1,3);
        queue(1,:) = [];
        neighbors = adj_list{node}.edges;
        imax = length(neighbors);
        for i = 1:imax
            if ~visited(neighbors(i))
                % Unvisited neighbor, add to queue and record it in output
                % variable. The edge used to reach this neighbor must be
                % unvisited too.
                comp_nodes{c}(end+1) = neighbors(i);
                visited(neighbors(i)) = true;
                % Bidirectional edge is visited now
                adj_visited{node}(i) = 1; 
                adj_visited{neighbors(i)}(adj_list{node}.rev_ind(i)) = 1;
                if ~isempty(trees{c})
                    % Create next tree depth info if needed
                    if length(trees{c}) < depth + 2
                        trees{c}{depth+2} = [];
                    end
                    % Add node info to the BFS tree
                    trees{c}{depth+2}(end+1).node = neighbors(i);
                    trees{c}{depth+2}(end).children = [];
                    trees{c}{depth+2}(end).par_ind = par_ind;
                    trees{c}{depth+2}(end).par_edge = i;
                    child_ind = length(trees{c}{depth+2});
                    trees{c}{depth+1}(par_ind).children(end+1) = child_ind;
                else
                    % The component is not a tree, this value does not
                    % matter anymore
                    child_ind = -1;
                end
                % Add neighbor to queue
                queue(end+1,:) = [neighbors(i), depth + 1, ...
                    child_ind]; %#ok<AGROW>
            elseif ~isempty(trees{c}) && adj_visited{node}(i) == 0
                % Unvisited edge to goes to a visited node: we have a cycle
                trees{c} = [];
            end
        end
    end
    
    % If the source is isolated, remove it from output vector
    if length(comp_nodes{c}) <= 1
        comp_nodes(c) = [];
        trees(c) = [];
        c = c - 1;
    end
    
    % BFS finished, find next unvisited node. If it exists then it must be
    % in a different component.
    while source_ind <= N
        source = source_list(source_ind);
        if ~visited(source), break; end
        source_ind = source_ind + 1;
    end
end

