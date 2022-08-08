function [ep, sol] = calculate_ep(ep_ocp, x0_p, x0, w_ep, costfun, sol_init)
arguments
    ep_ocp (1,1) casadi.Opti
    x0_p casadi.MX
    x0 {mustBeNumeric}
    w_ep casadi.MX
    costfun casadi.MX
    sol_init (1,2) = ones(1,2)
end

ep_ocp.set_value(x0_p, x0)

% calculate extreme points
ep_ocp.set_value(w_ep, [1e0, 1e-3]);
if isa(sol_init, 'casadi.OptiSol')
     ep_ocp.set_initial( [ep_ocp.x; ep_ocp.lam_g], sol_init(1).value( [ep_ocp.x; ep_ocp.lam_g] ) )
end
sol(1) = ep_ocp.solve();
ep(1,:) = ep_ocp.value(costfun);
ep_ocp.set_value(w_ep, [1e-3, 1e0]);
if isa(sol_init, 'casadi.OptiSol')
     ep_ocp.set_initial( [ep_ocp.x; ep_ocp.lam_g], sol_init(2).value( [ep_ocp.x; ep_ocp.lam_g] ) )
end
sol(2) = ep_ocp.solve();
ep(2,:) = ep_ocp.value(costfun);

end