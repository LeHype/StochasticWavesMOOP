%Multiobjective OCP: Initialize OCP and set cost function and disturbance
timehorizon =50;
timestep = 0.4;
[ocp,x,u,d] = initializeOCP(timehorizon,timestep);
%Test Gitflow again

%% Specify Pareto algorythm 
Algo = 'nbi';   %Options are awds or nbi as of now
%%
ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));


%Take Both Costs

costs = ([-x(6,end) x(7,end)]);

%% Here the ocp is scalarized and then the chosen pareto Algorythm is run
switch Algo
    case 'nbi'
        [p_params, ep, norm_costfun] = scalarize_moocp( ocp, costs, method="nbi" );
        [sol, ~] = nbi2d( ocp,ep, p_params, 15);
    case 'awds'
        [p_params, ep, norm_costfun] = scalarize_moocp( ocp, costs, method="wmm" );
        p_params= transpose(p_params);
        [sol, ~] = awds( ocp, norm_costfun, ep, p_params );
    otherwise
        warning('No valid Algorythm selected')
end
    


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
time = [0:timestep:timehorizon];
i = 15;
figure(i+10)
subplot(2,2,1)
plot(time,sol(i).value(x(2,:)))
title(['Angle of Model'])
subplot(2,2,2)
plot(time,sol(i).value(u))
title(['Input u'])
subplot(2,2,3)
plot(time,sol(i).value(d))
title(['Wave exitation'])
subplot(2,2,4)
plot(time,sol(i).value(x(6,:)),'g')
ylabel(['Energy harvested'])
title(['Energy harvested and damage done'])
hold on 
yyaxis right
plot(time,sol(i).value(x(7,:)),'r')
ylabel(['Damage done'])
hold on
yyaxis right

