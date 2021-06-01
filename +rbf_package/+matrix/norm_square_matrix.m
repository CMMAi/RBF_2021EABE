function [ m_distance ] = norm_square_matrix( measurement_points, basis_points )
%DISTANCE_MATRIX Summary of this function goes here
%   measurement_points : [M, dimension]
%   basis_points       : [N, dimension]
%   M: number of measurement points
%   N: number of basis points

[M, dimension] = size(measurement_points);
[N, dimension] = size(basis_points);

m_distance = zeros(M, N);

for dimension_index = 1:dimension
    [m_measurement, m_basis] = ndgrid...
        (measurement_points(:, dimension_index), ...
         basis_points(:, dimension_index));
     m_distance = m_distance + (m_measurement-m_basis).^2;
end

end

