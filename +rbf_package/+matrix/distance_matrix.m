function [ m_distance ] = distance_matrix( measurement_points, basis_points )
%DISTANCE_MATRIX Summary of this function goes here
%   measurement_points : [M, dimension]
%   basis_points       : [N, dimension]
%   M: number of measurement points
%   N: number of basis points

import rbf_package.matrix.norm_square_matrix;

m_distance = norm_square_matrix( measurement_points, basis_points );

m_distance = sqrt(m_distance);

end

