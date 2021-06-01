clear;

import rbf_package.*;
import rbf_package.basis_function.NMQ.*;

%% PDEs
u = @(x, y, z) exp(x+y+z);
du_dx = @(x, y, z) exp(x+y+z);
du_dy = @(x, y, z)  exp(x+y+z);
du_dz = @(x, y, z)  exp(x+y+z);
biharmonic_u = @(x, y, z) 9*exp(x+y+z);
scale = 0.25;

%% Loading collocation points
skip = 7;
Ni = 1500;

%% collocation points
boundary_data = load('boundary.mat');
Gamma1 = boundary_data.Gamma1;
interior_data = load('interior.mat');
Omega = interior_data.Omega;

Nb = length(Gamma1);
N = Ni + Nb*2;

normal_data = load('normal.mat');
normal_Gamma1 = normal_data.normal_Gamma1;

Gamma2 = Gamma1;
normal_Gamma2 = normal_Gamma1;

x_i = [Omega; Gamma1; Gamma2];

%% Test points
test_data = load('test.mat');
test_points = test_data.test_points;

center = [0  -0.2 0];

%% exact
u_exact = u(test_points(:, 1), test_points(:, 2), test_points(:, 3));

%% create ghost point
ghost_data = load('ghost.mat');
x_c = ghost_data.x_c;

for R = [2,2.5,3,3.5,4,4.5,5,5.5]    
   %% scale Ghost point
    x_j = x_c *R + center;

    %% Matrix
    dm_eval = matrix.distance_matrix(test_points, x_j);
    dm_Omega = matrix.distance_matrix(Omega, x_j);
    dm_Gamma1 = matrix.distance_matrix(Gamma1, x_j);
    dm_Gamma2 = matrix.distance_matrix(Gamma2, x_j);

    xm_Gamma2 = matrix.difference_matrix(Gamma2(:, 1), x_j(:, 1));
    ym_Gamma2 = matrix.difference_matrix(Gamma2(:, 2), x_j(:, 2));
    zm_Gamma2 = matrix.difference_matrix(Gamma2(:, 3), x_j(:, 3));

    rhs = [biharmonic_u(Omega(:, 1), Omega(:, 2), Omega(:, 3));...
        u(Gamma1(:,1), Gamma1(:,2), Gamma1(:,3));...
          du_dx(Gamma2(:, 1), Gamma2(:, 2), Gamma2(:, 3)) .* normal_Gamma2(:, 1)...
        + du_dy(Gamma2(:, 1), Gamma2(:, 2), Gamma2(:, 3)) .* normal_Gamma2(:, 2)...
        + du_dz(Gamma2(:, 1), Gamma2(:, 2), Gamma2(:, 3)) .* normal_Gamma2(:, 3)];
   
    %% Shape parameter (the same length as basis)
    c_assesment = shape_parameter.modified_Franke(N, 2*R);
    [c, c_min, c_max] = shape_parameter.variable(c_assesment, length(x_j));

    %% computation
    em = rbf(c, dm_eval);
    cm_Omega = biharmonic_rbf_3D(c, dm_Omega);
    cm_Gamma1 = rbf(c, dm_Gamma1);
    
    drbf_dx_Gamma2 = drbf_dx(c, dm_Gamma2, xm_Gamma2);
    drbf_dy_Gamma2 = drbf_dx(c, dm_Gamma2, ym_Gamma2);
    drbf_dz_Gamma2 = drbf_dx(c, dm_Gamma2, zm_Gamma2);
    
    cm_Gamma2 = drbf_dx_Gamma2 .* normal_Gamma2(:, 1) ... 
              + drbf_dy_Gamma2 .* normal_Gamma2(:, 2) ... 
              + drbf_dz_Gamma2 .* normal_Gamma2(:, 3);
    
    cm = [cm_Omega; cm_Gamma1; cm_Gamma2];
    
    %% solve
    u_solve = em*(cm\rhs);    
    diff = u_exact - u_solve;
    error = norm(diff,Inf); %diff.^2;

    fprintf('R =%3.1f, Max error = %10.3e, c_ini=%5.3f, c_range=(%5.3f, %5.3f)\n',...
        R, max(error), c_assesment, c_min, c_max );
end