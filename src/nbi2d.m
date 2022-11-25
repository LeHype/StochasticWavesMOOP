function [sol, p_param_vals] = nbi2d(ocp, ep, p_params, nPoints)
% nbi2d applies the normal-boundary intersection to an optimization problem:
%
% INPUTS:
% ocp: casadi.Opti object containing a scalarized optimization problem.
%
% ep: Extreme points of the optimization problem described by ocp.
%
% p_params: Cell array of the two Pareto parameter vectors, containing the
% parameters search direction and starting point.
%
% nPoints: Number of Pareto optimal points to be calculated.

arguments
    ocp (1,1) casadi.Opti
    ep  {mustBeNumeric}
    p_params cell
    nPoints (1,1) {mustBeNumeric}
end

% Weights are used in this context to describe the boundary plane as the
% convex combinations of all extreme points
weights = linspace(1,0,nPoints);
weights = [weights; 1-weights];
[spt, sdir] = p_params{:}; % get starting point and search direction parameters

% The search direction is the normal vector to the boundary plane
ndir = null(ep(2:end,:)-ep(1,:));
if all(ndir >= 0)
    ndir = -ndir;
end
ocp.set_value(sdir, ndir);

bp_pts = weights(1,:)'.*ep(1,:) + weights(2,:)'.*ep(2,:);

for iPt = 1:nPoints
    if iPt > 1
        % initialize the new ocp with a previous solution
        ocp.set_initial([ocp.x; ocp.lam_g], sol(iPt-1).value([ocp.x; ocp.lam_g]))
    end
    ocp.set_value( spt, bp_pts(iPt, :) )
    sol(iPt) = ocp.solve();
end

p_param_vals = bp_pts;

end

