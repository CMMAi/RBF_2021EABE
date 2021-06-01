function [ points ] = disk_2d( D, domain_shift, number_of_points )
%DISK_2D_POINT_GENERATOR Summary of this function goes here
%   Detailed explanation goes here
%   D                : domain diameter
%   domain_shift     : center shift of domain

N = 100; % number of boundary points
R = 0.5 * D; % domain radius;

t = 2 * pi * (1:N)/N;
domain = [R*cos(t') R*sin(t')];
domain = domain + domain_shift;
%plot(domain(:, 1), domain(:, 2), 'b'); hold on;
% create halton points
n_samples = 100000;
p_set = haltonset(2, 'Skip', floor(1000), 'Leap', 1e2);

samples = net(p_set, n_samples);
samples = samples*2*R - R + domain_shift;

is_within_domain = inpolygon(samples(:,1), samples(:,2), domain(:,1), domain(:,2));
samples = samples(is_within_domain, :);

points = samples(1:number_of_points, :);
end

