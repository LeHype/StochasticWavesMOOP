function [ocp, x, u, varargout] = ode2ocp(odefun, nx, nu, nSteps, timing, args)
% ode2ocp transforms a given ode with nx states and nu inputs to an ocp in
% an casadi.Opti() environment.
% INPUTS:
%   odefun: rhs of x' = f(x,u)
%
%   nx: number of states
%
%   nu: number of inputs
%
%   nSteps: Number of steps the system is predicted into the future
%
%   timing: Is either the temporal step size dt or 'optimal'. If 'optimal,
%   the time T will be a decision variable and be given out as the 4th
%   output.
%
% OPTIONAL ARGUMENTS:
%   x0: Numeric value for x0 to be fixed or 'param' to make it an Opti
%   parameter and make it an output.
%
%   xf: Numeric value for xf to be fixed or 'param' to make it an Opti
%   parameter and make it an output.
%
%   x_box: Box constraints on the states, has to be given in size (nx,2).
%   Nonexistent bounds can be given with -inf and inf.
%
%   u_box: Box constraints on the inputs, has to be given in size (nu,2).
%   Nonexistent bounds can be given with -inf and inf.
%
%   solver: NLP solver to be used, standard is 'ipopt'.
%
%   nd: Number of disturbances.
%
%   foh: Use first-order hold for inputs and disturbances
%
%   verbose: Verbosity level of the solver, currently only working for
%   ipopt.
%
% OUTPUTS:
%   ocp: casadi.Opti() object for the OCP
%
%   x: Matrix consisting of Opti parameters, doubles and Opti variables
%   based on the settings.
%
%   u: Opti.variables matrix for the inputs for the horizon.
%   
% OPTIONAL OUTPUTS: 
%   T: The length of the horizon as a decision variable and Opti variable.
%
%   x0: The parameter for the initial state x0, if x0 is 'param'.
%
%   xf: The parameter for the final state xf, if xf is 'param'.
%
%   d: The parameters for the disturbance if nd is set.
%
%   [h_integ, h_zero]: Matrix of h(x) and the corresponding equality
%   constraints, latter are used to get dual variables.
%
%   [g_x, g_x_constr]: Matrix of g(x) and the corresponding inequality
%   constraints.
%   
%   [g_u, g_u_constr]: Matrix of g(u) and the corresponding inequality
%   constraints.


arguments
    odefun  function_handle
    nx      (1,1) {mustBeNumeric}
    nu      (1,1) {mustBeNumeric}
    nSteps  (1,1) {mustBeNumeric}
    timing  {mustBeTiming}
    args.x0      {} = double.empty(nx, 0)
    args.xf      {mustBeValidState(args.xf, nx)} = double.empty(nx, 0)
    args.x_box   (:,2) {mustBeOfDim(args.x_box, nx), mustBeIncreasing(args.x_box)} = [-inf(nx, 1), inf(nx, 1)]
    args.u_box   (:,2) {mustBeOfDim(args.u_box, nu), mustBeIncreasing(args.u_box)} = [-inf(nu, 1), inf(nu, 1)]
    args.solver  (1,:) {mustBeText} = 'ipopt'
    args.nd      (1,1) {mustBeNumeric} = 0
    args.foh     (1,1) logical = false
    args.verbose (1,1) {mustBeNumeric} = 1
    args.ds      (1,:) {mustBeText} = 'central'     %derivative scheme
end

foh = double(args.foh);

varargout = cell(1,7);
ocp = casadi.Opti();
solver_opts = struct;


if isequal(args.solver, "ipopt")
    solver_opts.print_level = args.verbose;
end
ocp.solver(args.solver, struct, solver_opts)

if isnumeric(timing)
    dt = timing;
elseif isequal(timing, 'optimal')
    T  = ocp.variable(1);
    varargout{1} = T;
    dt = T/nSteps;
end

if isnumeric(args.x0)
    x0 = args.x0;
else
    x0 = ocp.parameter(nx,1);
    varargout{2} = x0;
end
x0_fixed = ~isempty(x0);

if isnumeric(args.xf)
    xf = args.xf;
else
    xf = ocp.parameter(nx,1);
    varargout{3} = xf;
end
xf_fixed = ~isempty(xf);

x = [x0, ocp.variable(nx, nSteps+1 - (x0_fixed + xf_fixed)), xf];
u = [ocp.parameter(nu, foh), ocp.variable(nu, nSteps)];
if strcmp(args.ds, 'backward')
du = [(u(2:end)-u(1:end-1))/dt (u(end)-u(end-1))/dt];
varargout{8} = du;
elseif strcmp(args.ds, 'forward')
du = [0 (u(2:end)-u(1:end-1))/dt];
varargout{8} = du;
elseif strcmp(args.ds, 'central')
du = [ 0 ((u(3:end)-u(2:end-1))/dt +  (u(2:end-1)-u(1:end-2))/dt)/2  ((u(end)-u(end-1))/dt +  (u(end-1)-u(end-2))/dt)/2 ];
varargout{8} = du;
else 
    warning('please choose suitable derivative method [backward ,forward, central]')
end


if args.nd > 0
    d = ocp.parameter(args.nd, nSteps + foh);
    varargout{4} = d;
end


h_integ = casadi.MX(nx, nSteps);
if args.nd > 0
    for iStep = 1:nSteps
        h_integ(:, iStep) = integrator_step_disturbed_du(x(:,iStep), u(:, iStep + (0:foh)), dt, odefun, d(:, iStep + (0:foh)), du(:, iStep + (0:foh))) - x(:,iStep+1);
        ocp.subject_to( )
    end
else
    for iStep = 1:nSteps
        h_integ(:, iStep) = integrator_step(x(:,iStep), u(:, iStep + (0:foh)), dt, odefun) - x(:,iStep+1);
    end
end
h_zero = h_integ(:) == 0;
varargout{5} = [h_integ(:), h_zero];
ocp.subject_to( h_zero )

if any(~isinf(args.x_box), 'all')
    is_inf = isinf(args.x_box);
    g_x_lower = args.x_box(~is_inf(:,1), 1) - x(find(~is_inf(:,1)), (1+x0_fixed):(end-xf_fixed));
    g_x_upper = x(find(~is_inf(:,2)), (1+x0_fixed):(end-xf_fixed)) - args.x_box(~is_inf(:,2), 2);
    g_x = [g_x_lower(:); g_x_upper(:)];
    g_x_constr = g_x <= 0;
    ocp.subject_to( g_x_constr )
    varargout{6} = [g_x, g_x_constr];
else
    varargout{6} = double.empty(0,2);
end

if any(~isinf(args.u_box), 'all')
    is_inf = isinf(args.u_box);
    g_u_lower = args.u_box(~is_inf(:,1), 1) - u(find(~is_inf(:,1)), :);
    g_u_upper = u(find(~is_inf(:,2)), :) - args.u_box(~is_inf(:,2), 2);
    g_u = [g_u_lower(:); g_u_upper(:)];
    g_u_constr = g_u <= 0;
    ocp.subject_to(g_u_constr)
    varargout{7} = [g_u, g_u_constr];
else
    varargout{7} = double.empty(0,2);
end

end

%%
function x_end = integrator_step(x0, u, dt, odefun)
% calculate one integration step with step size dt
import casadi.*

x0_rk = x0;
k = casadi.MX( size( x0, 1 ), 4 );

if size(u,2) == 1
    k(:,1) = odefun(x0_rk(:,end)                  , u);
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u);
else
    k(:,1) = odefun(x0_rk(:,end)                  , u(:,1));
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u*[0.5; 0.5]);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u*[0.5; 0.5]);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u(:,2));
end

x_end  = x0_rk + dt / 6 * k * [1 2 2 1]';

end

%%
function x_end = integrator_step_disturbed(x0, u, dt, odefun, d)
% calculate one integration step with step size dt
import casadi.*

x0_rk = x0;
k = casadi.MX( size( x0, 1 ), 4 );

if size(u,2) == 1
    k(:,1) = odefun(x0_rk(:,end)                  , u, d);
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u, d);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u, d);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u, d);
else
    k(:,1) = odefun(x0_rk(:,end)                  , u(:,1),       d(:,1));
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u(:,2),       d(:,2));
end

x_end  = x0_rk + dt / 6 * k * [1 2 2 1]';

end

%%
function x_end = integrator_step_disturbed_du(x0, u, dt, odefun, d, du)
% calculate one integration step with step size dt
import casadi.*

x0_rk = x0;
k = casadi.MX( size( x0, 1 ), 4 );

if size(u,2) == 1
    k(:,1) = odefun(x0_rk(:,end)                  , u, d,du);
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u, d,du);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u, d,du);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u, d,du);
else
    k(:,1) = odefun(x0_rk(:,end)                  , u(:,1),       d(:,1),       du(:,1));
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u*[0.5; 0.5], d*[0.5; 0.5], du*[0.5; 0.5]);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u*[0.5; 0.5], d*[0.5; 0.5], du*[0.5; 0.5]);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u(:,2),       d(:,2),       du(:,2));
end

x_end  = x0_rk + dt / 6 * k * [1 2 2 1]';

end

%%
function mustBeOfDim(input,dim)
    % Test for number of dimensions    
    if ~(size(input,1) == dim)
        eid = 'Size:wrongDimensions';
        msg = "Input must have size (" + dim + ",:).";
        throwAsCaller(MException(eid,msg))
    end
end

%%
function mustBeValidState(input, dim)
    if ~isnumeric(input) && ~isequal(input, 'param') || isnumeric(input) && ~(size(input,1) == dim)
        eid = 'Type:notState';
        msg = "Input must be either a positiv valid state of size (" + nx + ",1) or 'param'.";
        throwAsCaller(MException(eid,msg))
    end
end

%%
function mustBeTiming(input)
 % Test for valid timing
    if (~isnumeric(input) && ~isequal(input, 'optimal')) || (isnumeric(input) && input <= 0)
        eid = 'Type:wrongTiming';
        msg = "Input must be either a positiv scalar or 'optimal'.";
        throwAsCaller(MException(eid,msg))
    end
end

%%
function mustBeIncreasing(input)
 % Test for valid timing
    if ~all(input(:, 2:end) - input(:, 1:(end-1)) >= 0)
        eid = 'Value:decreasing';
        msg = "Input must be given in increasing order'.";
        throwAsCaller(MException(eid,msg))
    end
end