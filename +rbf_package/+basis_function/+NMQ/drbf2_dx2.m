function [ value ] = drbf2_dx2( c, r, x )
%DRBF2_DX2 Summary of this function goes here
%   Detailed explanation goes here
import rbf_package.basis_function.NMQ.rbf;

value = c.^2.*(rbf(c, r).^2 - c.^2.*x.^2)./rbf(c, r).^3;

end

