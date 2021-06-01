function [ value ] = biharmonic_rbf_3D( c, r )
%   This is NMQ RBF
import rbf_package.basis_function.NMQ.rbf;
value =  -(15*c.^4)./rbf(c, r).^7;
end

