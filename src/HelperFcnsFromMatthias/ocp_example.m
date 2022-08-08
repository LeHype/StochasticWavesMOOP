%% Building the MOOCP
% a priori known constraints
x0 = [0 0 0 225/180*pi 0 0]';
xf = [7, 5, 0, pi/2, 0]';
x_box = [-Inf 7.1; -Inf, Inf; -0.2 Inf; -Inf, Inf; -45/180*pi, 45/180*pi; -Inf, Inf];
u_box = [-10, 5; -20/180*pi, 20/180*pi];

% build ocp from ode using rk4, inputs are the ode, number of states nx,
% number of inputs nu, number of discretization steps nSteps, time
% discretization dt and some constraints
[ocp, x, u, T, ~, ~] = ode2ocp(@robot_ode, 6, 2, 60, 'optimal', x0=x0, x_box=x_box, u_box=u_box, foh=true);

% add additional constraints
ocp.subject_to( 1 <= T <= 20)
ocp.subject_to( x(1:5,end) == xf )
costfun = [T, x(6, end)];

%% Solve MOOCP
% scalarize the moocp to a parameterized siocp, possible methods are the
% weighted sum method ("ws"), weighted min max method ("wmm") and
% normal-boundary intersection/Pascoletti-Serafini ("nbi"/"ps").
[p_params, ep, norm_costfun] = scalarize_moocp( ocp, costfun, method="wmm" );

% apply an iterative Pareto front scheme for 2 dimensions to the problem,
% only works with ws and wmm.
[sol, ~] = awds( ocp, norm_costfun, ep, p_params );

%% Alternative Workflow
[p_params, ocp, ep, norm_costfun] = scalarize_moocp( ocp, costfun, method="nbi" );
[sol, ~] = nbi2d( ocp, ep, p_params, 20 );

%% Plotting the Pareto Front
J = [];
for i = 1:length(sol)
    J = [J; sol(i).value(costfun)];
end

plot(J(:,1), J(:,2), '.')

%%
function xdot = robot_ode(x, u)
    l = 1;
    xdot = [x(3)*cos( x(4) + x(5) ); x(3)*sin( x(4) + x(5) ); u(1); x(3) / l * sin( x(5) ); u(2); 0.5*u(1)^2];
end