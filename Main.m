
%% Here the standard serialized version is to check the runtime. 

%Multiobjective OCP: Initialize OCP and set cost function and disturbance
timehorizon =100;
timestep = 0.4;
nPoints = 8;
[ocp,x,u,d] = initializeOCP(timehorizon,timestep);

% opti.set_initial(sol1.value_variables());

ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));


costs= ([x(6,end) x(7,end)]);
[p_params, ep, norm_costfun] = scalarize_moocp( ocp, costs, method="nbi" ,simplify=false);
%%
tic;
[sol, ~] = nbi2d( ocp,ep, p_params, nPoints);
disp('serial took')
toc

%%
J = [];

for i = 1:length(sol)
     J = [J; [sol(i).value(x(6,end)) sol(i).value(x(7,end))]];
end

scatter(J(:,1),J(:,2));
