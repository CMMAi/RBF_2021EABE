function [ c, c_min, c_max ] = variable( c_ini, size, lower_bound)
%VARIABLE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    lower_bound = 0.3;
end
c_min = max(lower_bound, c_ini-0.5);
c_max = c_ini+0.5;
samples = haltonset(1, 'Skip', 500);
c = c_min + (c_max-c_min)*net(samples, size)';

end

