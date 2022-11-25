function [sol, costs] = awds(ocp, costfun, ep, p_params, ptDistance)
% awds applies the recursive formulation of the adaptive
arguments
   ocp        (1,1) casadi.Opti 
   costfun    (1,:) casadi.MX
   ep         {mustBeNumeric}
   p_params   (:,1) casadi.MX
   ptDistance (1,1) {mustBeNumeric} = 0.1
end

[sol, costs] = applyAWDS(ocp, costfun, ep, ep, p_params, ptDistance, []);

end

function [sol, costs] = applyAWDS(ocp, costfun, parents, ep, p_params, ptDistance, parentalSol)

newWeight     = parents\ones(size(parents,2),1);
newWeight     = newWeight'/norm(newWeight);

if all(newWeight <= 0)
    newWeight = -newWeight;
elseif any(newWeight < 0) || min(newWeight)/max(newWeight) <= 1e-6
    sol = [];
    costs = [];
    return
end

ocp.set_value( p_params, newWeight );
if ~isempty(parentalSol)
    ocp.set_initial( [ocp.x; ocp.lam_g], parentalSol.value([ocp.x; ocp.lam_g]) )

%     solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
%             'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
%             'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8,'print_level',1);
%     ocp.solver('ipopt', struct(), solver_opt)
end
sol = ocp.solve();
child = sol.value(costfun);

costs = child;

recombinations = nchoosek( 1:(size(parents,1)), size(parents, 2) );

for iParents = 1:size(recombinations, 2)

    dr = distance_ratio(parents(recombinations(iParents),:), child, ep);

    if dr <= ptDistance || any(vecnorm((parents-child)./vecnorm(parents,2,2),2,2) <= ptDistance/10)
        continue
    end

    [solNew, costsNew] = applyAWDS(ocp, costfun, [child; parents(recombinations(iParents),:)],...
        ep, p_params, ptDistance, sol);

    costs = [costs; costsNew];
    sol = [sol, solNew];
end


end

%%
function dr = distance_ratio(previousPoints, newPoint, EP)
% calculates the distance ratio between previously found solutions and the
% new solution. in the case of n=2 dimensions, distance ratio is equal to
% the euclidian distance in normalized (decision) space.
numParents = size(previousPoints, 1);
r = max(EP,[],1) - min(EP, [], 1);
dr = 0;

for p=1:numParents
    dr = dr + norm( (newPoint - previousPoints(p, :))./ r(1,:) );
end

dr = dr/numParents;

end
