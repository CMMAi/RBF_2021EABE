function [ value ] = laplacian_rbf_2D( c, r )
% Laplace rbf of NMQ
import rbf_package.basis_function.NMQ.rbf;
value =c.^2 .* (rbf(c, r) .^ 2 + 1) ./rbf(c, r).^3;
end

