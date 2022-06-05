%Multiobjective OCP: Initialize OCP and set cost function and disturbance
timehorizon =50;
timestep = 0.4;
nPoints = 15
[ocp,x,u,d] = initializeOCP(timehorizon,timestep);

% opti.set_initial(sol1.value_variables());

ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));



[p_params,up,np, ep, norm_costfun] = Jank_scalarize_moocp( ocp, costs, method="nbi" );


%%
[ocp bp_pts,spt,ndir] = NewJank( ocp,ep, p_params, 15);
ocp.set_value(spt,bp_pts(1,:));


%%
J = zeros(nPoints,2);
Jprime = zeros(nPoints,2);
UARRAY = zeros((timehorizon/timestep)+1,nPoints);
parfor i = 1:15
    disp(i)
[ocp,x,u,d] = initializeOCP(timehorizon,timestep);


ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));

costs= ([-x(6,end) x(7,end)]-[up])./[np-up];


sdir = (ocp.parameter(2));
ocp.set_value(sdir, ndir);

l   = ocp.variable( 1 );
pt  = ocp.parameter( 1, 2 );
dir = ocp.parameter( 1, 2 );
ocp.subject_to( pt + l*dir >= costs )
        
ocp.minimize(-l)
ocp.set_value(pt,bp_pts(i,:));
ocp.set_value(dir,ndir);
ocp.solve()
J(i,:) = (ocp.value(costs));
Jprime(i,:) = ocp.value([x(6,end) x(7,end)]);
UARRAY(:,i) = ocp.value(u);
end

%%



tic;

timehorizon =50;
timestep = 0.4;
weight = [0:1/20:1];
solA = cell(20,1);
ocpA = cell(20,1);
UARRAY = zeros((timehorizon/timestep)+1,20);
for i = 1:20
    ocpA{i}=ocp.copy();
end
J = zeros(2,20);
parfor i = 1:20

 [ocpA{i},x,u,d] = initializeOCP(timehorizon,timestep);
 costs = [x(6,end) x(7,end)];
 ocpA{i}.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));
 ocpA{i}.minimize(costs(1)*weight(i)+costs(2)*(1-weight(i)));
 solA{i}=ocpA{i}.solve();

 UARRAY(:,i) = (solA{i}.value(u));
ocpA{i} = [];
 J(:,i) = [solA{i}.value(x(6,end)) solA{i}.value(x(7,end))];

end
disp('Computation took')
disp(toc);

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
% J = [];
J2 = [];
% for i = 1:length(sol)
%     J = [J; sol(i).value(costs)];
% end
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

