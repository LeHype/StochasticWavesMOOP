%Multiobjective OCP: Initialize OCP and set cost function and disturbance
timehorizon =20;
timestep = 0.1;
[ocp,x,u,d] = initializeOCP(timehorizon,timestep);

%%
ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));


%Take Both Costs
costs = ([-x(6,end) x(7,end)]);

%Scalarize i.e. get a sense of how much bigger one costfunction is then the
%other
[p_params, ep, norm_costfun] = scalarize_moocp( ocp, costs, method="wmm" );
p_params= transpose(p_params);
%% Run the multiobjective obtimization
[sol, ~] = awds( ocp, norm_costfun, ep, p_params );
%% Plotting the Pareto Front  \\This is super jank 
J = [];
J2 = [];
for i = 1:length(sol)
    J = [J; sol(i).value(costs)];
end
NonDomInc = ParetoFilter(J);
for i = 1:length(NonDomInc)
    J2(i,:) = [[J(NonDomInc(i),:)],NonDomInc(i)];
end
figure
for i = 1:length(J2)
    scatter(J2(i,1), J2(i,2), 100,'filled')
    text(J2(i,1), J2(i,2),['\leftarrow',int2str(J2(i,3))],'Interpreter','tex')
    hold on 
end

%%



%% Plot if wanted
i = 16;
figure(i+10)
subplot(2,2,1)
plot(sol(i).value(x(2,:)))
title(['Angle of Model'])
subplot(2,2,2)
plot(sol(i).value(u))
title(['Input u'])
subplot(2,2,3)
plot(sol(i).value(d))
title(['Wave exitation'])
subplot(2,2,4)
plot(sol(i).value(x(6,:)),'g')
ylabel(['Energy harvested'])
title(['Energy harvested and damage done'])
hold on 
yyaxis right
plot(sol(i).value(x(7,:)),'r')
ylabel(['Damage done'])
hold on
yyaxis right

