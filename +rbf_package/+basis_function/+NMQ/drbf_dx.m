function [ value ] = drbf_dx( c, r, x )
%DRBF_DX Summary of this function goes here
%   Detailed explanation goes here
import rbf_package.basis_function.NMQ.rbf;
value = c.^2.*x./rbf(c, r);

end

