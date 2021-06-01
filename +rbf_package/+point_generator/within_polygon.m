function [ points, D, center ] = within_polygon( polygon_points, number_of_points )
%WITHIN_POLYGON_POINT_GENERATOR Summary of this function goes here
%   Detailed explanation goes here
% polygon_points: polygon described with points in [x1, y1; x2, y2;....; xn, yn]

import rbf_package.point_generator.disk_2d;

center = [mean(polygon_points(:, 1)), mean(polygon_points(:, 2))];

x_d = abs(max(polygon_points(:, 1)) - min(polygon_points(:, 1)));
y_d = abs(max(polygon_points(:, 2)) - min(polygon_points(:, 2)));

D = norm([x_d, y_d]);
R = 0.5*D;

%% create circle domain
samples = disk_2d(D, center, 10000);


%% generate points
is_within_domain = inpolygon(samples(:,1), samples(:,2), polygon_points(:,1), polygon_points(:,2));

samples = samples(is_within_domain, :);

points = samples(1:number_of_points, :);

end

