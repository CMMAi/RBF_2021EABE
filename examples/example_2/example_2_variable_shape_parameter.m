clear;
import rbf_package.*;
import rbf_package.basis_function.NMQ.*;

%% modified equations
cm_form = @(c, r, xm, ym, xi, yi) laplacian_rbf_2D(c, r) + (xi.^2 + yi.^2) .* rbf(c, r)...
                        + yi .* cos(yi) .* drbf_dx(c, r, xm) + sinh(xi) .* drbf_dx(c, r, ym);
                    
%% PDEs
u = @(x, y) sin(pi*x).*cosh(y) - cos(pi*x).*sinh(y);
laplace_u = @(x, y) (cos(pi*x).*sinh(y) - sin(pi*x).*cosh(y))*(pi^2 - 1);
du_dx = @(x, y) pi*cos(pi*x).*cosh(y) + pi*sin(pi*x).*sinh(y);
du_dy = @(x, y) sin(pi*x).*sinh(y) - cos(pi*x).*cosh(y);
f = @(x, y) laplace_u(x, y) + (x.^2+y.^2) .* u(x, y) + ...
            y .* cos(y) .* du_dx(x, y) + sinh(x) .* du_dy(x, y);
gradient_u = @(x, y) [du_dx(x,y), du_dy(x,y)];

%% define number of points
Nb1 = 200;
Nb2 = Nb1;
Nb = Nb1 + Nb2;
Ni = 500;
N = Ni + Nb;

%% create collocation points
theta=linspace(0,2*pi,Nb+1); theta(end)=[];
theta1=theta(1:Nb1); theta2=theta(Nb1+1:end);

r_parametric = @(theta, scale) scale * ( 1+cos( 4 * theta).^2 );
r = @(theta) r_parametric(theta, 0.5);

Gamma1 = [r(theta1).*cos(theta1); r(theta1).*sin(theta1)]';
Gamma2 = [r(theta2).*cos(theta2); r(theta2).*sin(theta2)]';
Gamma = [Gamma1; Gamma2];

dx=gradient(Gamma(:,1));
dy=gradient(Gamma(:,2));
xnn=dy./sqrt(dx.^2+dy.^2);
ynn=-dx./sqrt(dx.^2+dy.^2); 
normal_Gamma2=[xnn(Nb1+1:end), ynn(Nb1+1:end)];

[Omega, ~, domain_center] = point_generator.within_polygon([Gamma1; Gamma2], Ni);

%% test points
t = linspace(0, 2*pi, 50);
gamma = [r(t).*cos(t); r(t).*sin(t)]';
[test, ~, ~] = point_generator.within_polygon(gamma, 100); %test points
test_point=[test;gamma];

%% exact
u_exact = u(test_point(:, 1), test_point(:, 2));
gradient_Gamma2 = gradient_u(Gamma2(:, 1), Gamma2(:, 2));
rhs = [f(Omega(:, 1), Omega(:, 2));...
        u(Gamma1(:,1), Gamma1(:,2));...
        gradient_Gamma2(:, 1).*normal_Gamma2(:, 1)+gradient_Gamma2(:, 2).*normal_Gamma2(:, 2)];

%% create ghost points
[xc, yc] = point_generator.fabric_pattern(N,1);
 
for R = [ 2, 2.5,3,3.5,4,4.5,5,5.5]
    %% diameter
    D = 2*R;
	%% scaling ghost points
	x_j = [xc*R yc*R];   

	%% Shape parameter (the same length as basis)
    c_assesment = shape_parameter.modified_Franke(N, 2*R);
    [c, c_min, c_max] = shape_parameter.variable(c_assesment, length(x_j));

	%% Matrix
	dm_eval = matrix.distance_matrix(test_point,x_j);
    dm_Omega = matrix.distance_matrix(Omega, x_j);
    dm_Gamma1 = matrix.distance_matrix(Gamma1, x_j);
    dm_Gamma2 = matrix.distance_matrix(Gamma2, x_j);

    xm_Omega = matrix.difference_matrix(Omega(:, 1), x_j(:, 1));
    ym_Omega = matrix.difference_matrix(Omega(:, 2), x_j(:, 2));

    xm_Gamma1 = matrix.difference_matrix(Gamma1(:, 1), x_j(:, 1));
    ym_Gamma1 = matrix.difference_matrix(Gamma1(:, 2), x_j(:, 2));

    xm_Gamma2 = matrix.difference_matrix(Gamma2(:, 1), x_j(:, 1));
    ym_Gamma2 = matrix.difference_matrix(Gamma2(:, 2), x_j(:, 2));

    %% Assembly
    em = rbf(c, dm_eval);
    cm_Omega = cm_form(c, dm_Omega, xm_Omega, ym_Omega, Omega(:, 1), Omega(:, 2));
    cm_Gamma1 = rbf(c, dm_Gamma1);

    drbf_dx_Gamma2 = drbf_dx(c, dm_Gamma2, xm_Gamma2);
    drbf_dy_Gamma2 = drbf_dx(c, dm_Gamma2, ym_Gamma2);
    cm_Gamma2 = drbf_dx_Gamma2 .* normal_Gamma2(:, 1) + drbf_dy_Gamma2 .* normal_Gamma2(:, 2);

    cm = [cm_Omega; cm_Gamma1; cm_Gamma2];

    %% solve
    coefficients = cm\rhs;
    u_solve = em*coefficients;    
    diff = u_exact - u_solve;
    error = norm(diff, Inf);

    fprintf('R =%3.1f,  Maxerr = %10.3e, c_i=%5.3f, c_range=(%5.3f, %5.3f)\n',...
	R, error, c_assesment, c_min, c_max);

end
