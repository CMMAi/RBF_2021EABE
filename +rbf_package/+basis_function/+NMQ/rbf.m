function [ value ] = rbf( c, r )
%   This is NMQ RBF

value = sqrt(1+(c.*r).^2);

end

