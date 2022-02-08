function p = radix_sort(v)
%RADIX_SORT Radix sort for unsigned 64-bit integers.
%   p = pigeonhole(v)
%   Sorts the L entries of vector v in increasing order. The entries of v
%   are assumed to be integers. 
%   RADIX_SORT returns a permutation vector p such that v(p) = v_sorted. 
%   This sorting algorithm is stable and runs in O(L) time for unsigned 
%   64-bit integers.


%% Write integer in radix base (radix=256), O(N)
% First we remove negative entries
v = uint64(v);
v_min = min(v);
if v_min < 0
    v = v - v_min;
end

% Digit matrix
dv = zeros(length(v), 8, 'uint16');
L = length(v);
for i = 1:8
    dv(:,i) = uint16(bitshift(bitshift(v, 8.*(i-1), 'uint64'), ...
        -56, 'uint64'));
end


%% Counting sort subroutine, O(N)
p = (1:L).';
for i = 8:-1:1
    p_i = counting_sort(dv(p,i)+1, 256);
    p = p(p_i);
end

