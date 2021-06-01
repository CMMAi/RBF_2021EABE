% Poison equation
% comparason between ghost point and Kansa method
clear;
warning off all;
import rbf_package.*;
import rbf_package.basis_function.NMQ.*;

u = @(x, y) exp(x + y);
laplace_u = @(x, y) 2*exp(x + y);

%% define number of points
Nb = 100;
N = 400;
Ni = N - Nb;

theta = linspace(0, 2*pi, Nb+1);
theta(end)=[];
r_parametric = @(theta, scale) scale.*(cos(3*theta) + (2 - sin(3*theta).^2).^0.5).^(1/3);
r = @(theta) r_parametric(theta, 1);

Gamma = [r(theta).*cos(theta); r(theta).*sin(theta)]';
[Omega, ~, domain_center] = point_generator.within_polygon(Gamma, Ni);
x_i = [Omega; Gamma];


%% test points
t = linspace(0, 2*pi, 50);
gamma = [r(t).*cos(t); r(t).*sin(t)]';
[test, ~, ~] = point_generator.within_polygon(gamma, 100); %test points
test_point=[test;gamma];

%% exact
u_exact = u(test_point(:, 1), test_point(:, 2));
rhs = [laplace_u(Omega(:, 1),Omega(:, 2)); u(Gamma(:,1), Gamma(:,2))];

%% ghost points
[xc, yc] = point_generator.fabric_pattern(N,1);

for R = [ 2, 2.5,3,3.5,4,4.5,5,5.5]
   
	%% scale ghost points
    x_j =[R*xc, R*yc];

    %% Matrix
    dm_eval = matrix.distance_matrix(test_point, x_j);
    dm_Omega = matrix.distance_matrix(Omega, x_j);
    dm_Gamma = matrix.distance_matrix(Gamma, x_j);

    %% Shape parameter
    c = shape_parameter.modified_Franke(N, 2*R);
    
    %% computation
    em = rbf(c, dm_eval);
    cm_Omega = laplacian_rbf_2D(c, dm_Omega);
    cm_Gamma = rbf(c, dm_Gamma);
    cm = [cm_Omega; cm_Gamma];    
    %% solve
    u_solve = em*(cm\rhs);    
    diff = u_exact - u_solve;
    error = norm(diff,Inf);
    
    fprintf('R =%4.1f,  Maxerr = %10.3e,  c=%5.3f\n', R, max(error), c);
end
