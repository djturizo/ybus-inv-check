function p = counting_sort(v, N)
%COUNTING_SORT Counting sort for vectors of integer entries.
%   p = counting_sort(v, N)
%   Sorts the L entries of vector v in increasing order. The entries of v
%   are assumed to be integers in the range [1,N]. 
%   COUNTING_SORT returns a permutation vector p such that v(p) = v_sorted. 
%   This sorting algorithm is stable and runs in O(L+N) time.


%% Initialize variables, O(L+N)
count = zeros(N, 1);
L = length(v);
p = zeros(L, 1);


%% Main loop, O(L+N)
% Histogram count
for k = 1:L
    count(v(k)) = count(v(k)) + 1;
end

% Prefix sum
count = cumsum(count);

% Build permutation vector
for k = L:-1:1
    count(v(k)) = count(v(k)) - 1;
    p(count(v(k))+1) = k;
end

