function [ value ] = biharmonic_rbf( c, r )
%   This is NMQ RBF
import rbf_package.basis_function.NMQ.rbf;
value = (c.^4.*(c.^4.*r.^4 + 8*c.^2.*r.^2 - 8))./(c.^2.*r.^2 + 1).^(7/2);

end

