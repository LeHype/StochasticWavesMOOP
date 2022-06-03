%Base Version: Initialize OCP and set cost function and disturbance
timehorizon =100;
timestep = 0.1;
[ocp,x,u,d] = initializeOCP(timehorizon,timestep);

%%
ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));
ocp.minimize(x(6,end));
sol = ocp.solve();

%% Plot if wanted

figure(1)
plot(sol.value(x(2,:)))
hold on
yyaxis right
plot(sol.value(u))