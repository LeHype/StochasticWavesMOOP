x_val  = [1 1 0 0]';
u_val = 0;
xf = [0,0]';
nHorizon = 60;
nPoints = 20;
foh = true;
normalize = "fix";
% weight = [0.0514 0.9987];
weight = [0.5 0.5];

u_box = [-1, 1];

[ocp, x, u, ~, x0_p, ~] = ode2ocp(@oszi_ode, 4, 1, nHorizon, 0.1, x0='param', u_box=u_box, foh=foh);

accum_costs = zeros(1,2);
costfun = [x(3,end), x(4,end)];

ocp.set_value(x0_p, x_val(:,end))
if foh
    ocp.set_value(u(:,1), u_val(:,end))
end
if isequal(normalize, "adapt")
    [p_params, up, np, norm_costfun, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize=normalize );
elseif isequal(normalize, "none")
    [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize=normalize );
elseif isequal(normalize, "fix")
    [p_params, ep_val, norm_costfun, ep, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize=normalize );
end

sol = {};
J = {};
sol_init = [];
idx = uint8.empty(0,nHorizon);

for i = 1:nHorizon
    ocp.set_value(x0_p, x_val(:,end))
    if foh
        ocp.set_value(u(:,1), u_val(:,end))
    end
    if i > 1
        ocp.set_initial(x(:,2:end-1), sol{i-1}(idx(i-1)).value(x(:,3:end)));
        ocp.set_initial(u(:,1+foh:end-1), sol{i-1}(idx(i-1)).value(u(:,2+foh:end)));
    end

    if ismember(normalize, ["adapt", "none"])
        ep_fun_in = {ep_ocp, x0_p, x_val(:, end), w_ep, costfun};
        if i > 1
            ep_fun_in{end+1} = sol_ep;
        end
        [ep_val, sol_ep] = calculate_ep(ep_fun_in{:});
    end
    
    if isequal(normalize, "fix")
        up_val = min(ep);
        np_val = max(ep);
    elseif isequal(normalize, "none")
        up_val = zeros(1,2);
        np_val = ones(1,2);
    else
        up_val = min(ep_val);
        np_val = max(ep_val);
    end

    if isequal(normalize, "adapt")
        ocp.set_value(up, up_val);
        ocp.set_value(np, np_val);
        ep_val_norm = rescale_upnp(ep_val, up_val, np_val);
    elseif isequal(normalize, "none")
        norm_costfun = costfun;
        ep_val_norm = ep_val;
    else
        ep_val_norm = rescale_upnp(ep_val, up_val, np_val);
    end

    [sol{i}, ~] = nbi2d( ocp, ep_val_norm, p_params, nPoints );
    J_cell = arrayfun(@(s) s.value(costfun), sol{i}, 'uni', 0);
    J{i} = vertcat(J_cell{:});
    [~,idx(i)] = min( rescale_upnp(J{i}, up_val, np_val)*weight' );

    sol_init = sol{i}(idx(i));
    x_val(:, end+1) = sol{i}(idx(i)).value(x(:,2));
    u_val(:, end+1) = sol{i}(idx(i)).value(u(:,1+foh));

    ep_ocp.set_initial(x(:,2:end), sol_init.value(x(:,[3:end, end])))
    ep_ocp.set_initial(u(:,1+foh:end), sol_init.value(u(:,[2+foh:end, end])))
    plot(x_val(1,:), x_val(2,:), sol{i}(idx(i)).value(x(1,2:end)), sol{i}(idx(i)).value(x(2,2:end)), "--", 'LineWidth', 2)
end

%%
function xdot = oszi_ode(x, u)
    xdot = [x(2); -x(1) + u; x(1:2)'*x(1:2); 0.5*u^2];
end

%%
function out = rescale_upnp(vec, up, np)

out = (vec - up)./(np - up);

end