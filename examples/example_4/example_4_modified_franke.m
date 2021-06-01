clear;
import rbf_package.*;
import rbf_package.basis_function.NMQ.*;

%% Material
q = -1e6;
E = 2.1e11;
h = 0.01;
poisson_ratio = 0.3;
r0 = 1;

material_tangent = (E*h^3)/(12*(1-poisson_ratio^2));

force = q/material_tangent;

%% exact
to_x = @(r, theta) r .* cos(theta);
to_y = @(r, theta) r .* sin(theta);

to_polar = @(x, y) [sqrt(x.^2 + y.^2),  atan(y/x)];
to_polar_r = @(x, y) sqrt(x.^2 + y.^2);
to_polar_theta = @(x, y) atan(y/x);

u = @(x, y) (1./64) .* (q/material_tangent) .* (r0^2 - to_polar_r(x, y).^2).^2;

%% define number of points
Nb = 100;
Ni = 300;
N = Ni + Nb + Nb;

%% Collocation points
theta = linspace(0, 2*pi-1e-10, Nb);
center = [0, 0];

r_parametric = @(theta, scale) r0;
r = @(theta) r_parametric(theta, 1);

Gamma1 = [r(theta).*cos(theta); r(theta).*sin(theta)]';

gradient_Gamma1 = [gradient(Gamma1(:, 1)), gradient(Gamma1(:, 2))];
normal_Gamma1 = [gradient_Gamma1(:, 2), -gradient_Gamma1(:, 1)] ./ sqrt(gradient_Gamma1(:, 1).^2 + gradient_Gamma1(:, 2).^2);

Gamma2 = Gamma1;
gradient_Gamma2 = gradient_Gamma1;
normal_Gamma2 = normal_Gamma1;

[Omega, ~, domain_center] = point_generator.within_polygon(Gamma1, Ni);

x_i = [Omega; Gamma1; Gamma2];

%% test points
t = linspace(0, 2*pi, 133);
gamma = [r(t).*cos(t); r(t).*sin(t)]';
[test, ~, ~] = point_generator.within_polygon(gamma, 400); %test points
test_point=[test;gamma]; %nt=length(test(:,1));

%% exact
u_exact = u(test_point(:, 1), test_point(:, 2));

%% generate ghost point
[xc, yc] = point_generator.fabric_pattern(N,1);

for  R = [2,2.5,3,3.5,4,4.5,5,5.5]
       
    %% scale ghost point
    x_j = [xc*R yc*R];

    %% Matrix   
    dm_eval = matrix.distance_matrix(test_point, x_j);
    dm_Omega = matrix.distance_matrix(Omega, x_j);
    dm_Gamma1 = matrix.distance_matrix(Gamma1, x_j);
    dm_Gamma2 = matrix.distance_matrix(Gamma2, x_j);

    xnm_Gamma2 = repmat(normal_Gamma2(:, 1), 1, N);
    ynm_Gamma2 = repmat(normal_Gamma2(:, 2), 1, N);

    xm_Gamma2 = matrix.difference_matrix(Gamma2(:, 1), x_j(:, 1));
    ym_Gamma2 = matrix.difference_matrix(Gamma2(:, 2), x_j(:, 2));

    rhs = [ones(Ni, 1) .* force ; zeros(Nb, 1); zeros(Nb, 1)];

    %% Shape parameter (the same length as basis)
    c = 0.8 * N^.25/(2*R);
    
    %% Assembly
    em = rbf(c, dm_eval);
    cm_Omega = biharmonic_rbf_2D(c, dm_Omega);
    cm_Gamma1 = rbf(c, dm_Gamma1);
        
    drbf_dx_Gamma2 = drbf_dx(c, dm_Gamma2, xm_Gamma2);
    drbf_dy_Gamma2 = drbf_dx(c, dm_Gamma2, ym_Gamma2);
    cm_Gamma2 = drbf_dx_Gamma2 .* xnm_Gamma2 + drbf_dy_Gamma2 .* ynm_Gamma2;
    
    cm = [cm_Omega; cm_Gamma1; cm_Gamma2];
    
    %% solve
    u_solve = em*(cm\rhs);    
    diff = u_exact - u_solve;
    error = norm(diff, Inf);
    
    fprintf('R =%3.1f, Maxerr = %10.3e, c=%5.3f\n', R, error, c);    
end