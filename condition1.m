function [cond_true, fail_node] = condition1(tree, adj_list, shunts, tol)
%CONDITION1 Evaluates Condition 1) for a tree power system
%   is_true = condition1(tree, adj_list, shunts, change_root)
%   Checks Condition 1) of Theorem 3 of the paper in a power system whose
%   network graph is desbribed by tree, and whose edges are a subset of
%   adj_list. tree must follow the structure of the output of the function
%   GET_COMPS, and adj_list must follow the structure of the output of 
%   BI_ADJ_LIST. adj_list may have more edges than the ones of tree. shunts
%   is a vector with the value of the shunt admittances of every node,
%   possibly including nodes not connected to tree. 
%   CONDITION1 returns the logic variable  is_true with a value of true if 
%   the condition holds or false otherwise.
%   CONDITION1 also returns fail_node, which is empty is the condition 
%   holds. Otherwise fail_node stores the value of the node where we found 
%   a violation.
%   If tree has N nodes and L edges, this algorithm runs in O(N+L) time.


%% Compute equivalent admittance going from leaves to root, O(N+L)
cond_true = true;
fail_node = [];
imax = length(tree);
for i = imax:-1:1
    jmax = length(tree{i});
    for j = 1:jmax
        % Compute y^sb_n fro Eq. (79) of the paper.
        ysb = shunts(tree{i}(j).node);
        hmax = length(tree{i}(j).children);
        for h = 1:hmax
            child = tree{i}(j).children(h);
            ysb = ysb + tree{i+1}(child).yb;
        end
        if i > 1
            % Get Y^p,n afrom Eq. (76) of the paper.
            par_ind = tree{i}(j).par_ind;
            par_edge = tree{i}(j).par_edge;
            parent = tree{i-1}(par_ind).node;
            Ypn = adj_list{parent}.Ybus{par_edge};
            if abs(Ypn(2,2) + ysb) <= tol
                % Condition 1) is violated at this node.
                cond_true = false;
            end
            % Compute y^b_n from Eq. (86) of the paper.
            tree{i}(j).yb = Ypn(1,1) - ...
                ((Ypn(1,2).*Ypn(2,1)) ./ (Ypn(2,2) + ysb));
        else
            % We are at the root, the equivalent admittance is just ysb.
            if abs(ysb) <= tol
                % Condition 1) violated at the root.
                cond_true = false;
            end
        end
        if ~cond_true
            fail_node = tree{i}(j).node;
            return;
        end
    end
end








