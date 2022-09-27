%% MCP WITH ClOSE LOOK AT THE DIFFERENT SOLUTION
%% START WITH NORMAL NBI TO FIND GOOD PARAMETERS FOR COST FUNCTION
timehorizon = 30;
timestep = 0.2; 
nPoints = 12;           %% Should be a fast algorythm so you can have lots of points
nCores = 6;             %% Please pay close attention to RAM. If full reduce 

solA = ParalellNBI(timehorizon,timestep,nPoints,nCores) ;

% delete(gcp('nocreate')) %% shutdown pool 

%% Plot results for convinient choice of PP
UU = zeros(length(solA),length(solA(1).u));
XX = zeros(length(solA),length(solA(1).x(2,:)));
DD = zeros(length(solA),length(solA(1).x(6,:)));
EE = zeros(length(solA),length(solA(1).x(6,:)));
PP= zeros(length(solA),2);
figure(1)

for i = 1:length(solA)
UU(i,:)=solA(i).u;
XX(i,:)=solA(i).x(2,:);
DD(i,:)=solA(i).x(7,:);
EE(i,:)=solA(i).x(6,:);
PP(i,:)=[solA(i).x(6,end) solA(i).x(7,end)];
% disp([num2str(solA(i).time),solA(i).warmstart])
scatter(PP(i,1),PP(i,2),100,'filled');
hold on
 text(PP(i,1), PP(i,2),['\leftarrow',int2str(i)],'Interpreter','tex');
end
title('Pareto Front')
xlabel(['Energy'])
ylabel('damage')

%% Further inquiries into the solution. Not part of the main workflow
warm = []
cold =[]
for i =1:length(solA)
disp([num2str(solA(i).time),'    ', solA(i).warmstart ])
if strcmp(solA(i).warmstart,'true')
    warm = [warm solA(i).time];
else
    cold = [cold solA(i).time];
end
end
mean(warm)
mean(cold)
%% Please chose a point from the pareto Front and then inspect it. 
Pindx = 1; %% Manually change this!
time = [0:timestep:timehorizon];
if length(time) < length(UU(Pindx,:));
    time = [time  time(end)]
end
figure(12)
subplot(2,2,1)
plot(time,UU(Pindx,:))
xlabel('time')
ylabel('Voltage')

subplot(2,2,2)
plot(time,XX(Pindx,:))
xlabel('time')
ylabel('angle')


subplot(2,2,3)
plot(time,XX(Pindx,:))
hold on
ylabel('angle')
yyaxis right
plot(time,UU(Pindx,:))
xlabel('time')
ylabel('voltage')

subplot(2,2,4)
plot(time,DD(Pindx,:))
hold on
ylabel('Damage')
yyaxis right
plot(time,EE(Pindx,:))
xlabel('time')
ylabel('Energy')
%% Start MCP here


[ocp,x,u,d,x0] = initializeOCP_NBI(timehorizon,timestep,solA(Pindx).ParetoParameters);

sol=ocp.solve();
uOld1= sol.value(u);
xOld1= sol.value(x);
dOld1= sol.value(d);
%% Start MPC
tic;
timeshift = 1;
timeshiftindx= max(1,round(length(x)/(timehorizon/timeshift)))+1;
x0new =xOld1(:,timeshiftindx);
x0new(6:7)= [0;0];
timeshift1 = time(timeshiftindx);
initialGuessX = [xOld1(:,timeshiftindx:end) repmat(xOld1(:,end),1,timeshiftindx-1)];
initialGuessU = [uOld1(timeshiftindx:end) repmat(uOld1(end),1,timeshiftindx-1)];

initialGuessX(6:7,:)=zeros(2,length(initialGuessX(1,:)));
xOld = xOld1;
uOld = uOld1;

timeold=time;
for i = 1:50





timeshift = timeshift1*i;
[ocp,x,u,d,x0] = initializeOCP_NBI(timehorizon,timestep,solA(Pindx).ParetoParameters,x0=x0new,ts=timeshift);
% ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep+timeshift:timestep:(d.length()*timestep)+timeshift]));
ocp.set_initial(x,initialGuessX);
ocp.set_initial(u(2:end),initialGuessU(2:end));
x0new(6:7) = [0 ; 0];
ocp.set_value(x(:,1),x0new)
ocp.set_value(u(1),uOld(timeshiftindx));
costs= ([x(6,end) x(7,end)]);
ocp.minimize( costs*[1 ; 1E-4]);
sol1 = ocp.solve();
ep(1,:) = sol1.value(costs);

ocp.minimize( costs*[1E-4 ; 1]);
sol2= ocp.solve();
ep(2,:) = sol2.value(costs);



up = min( ep, [], 1 );
np = max( ep, [], 1 );




costs= ([x(6,end) x(7,end)]-[up])./[np-up];


sol1 = ocp.solve();
time1 =[timestep+timeshift:timestep:(d.length()*timestep)+timeshift];
f = figure(1);
subplot(2,1,1)
plot(timeold,uOld,"Color",'r','LineStyle','-','LineWidth',5)
hold on 
 plot(time1,sol1.value(u),"Color",'b','LineStyle','--','LineWidth',5)
hold off
ylim([0 1200]);
xlim([timeold(1) time1(end)]);
subplot(2,1,2)
plot(timeold,xOld(2,:),"Color",'r','LineStyle','-','LineWidth',5)
hold on 
plot(time1,sol1.value(x(2,:)),"Color",'b','LineStyle','--','LineWidth',5)
hold off
ylim([-1.5 1.5]);
xlim([timeold(1) time1(end)]);
saveU = getframe(f).cdata;
cd Images/
filename =[num2str(i),'_U.png'];
imwrite(saveU ,filename);
cd ..
% ff = figure(2);
% 
% plot(timeold,xOld(2,:),"Color",'r','LineStyle','-','LineWidth',5)
% hold on 
% plot(time1,sol1.value(x(2,:)),"Color",'b','LineStyle','--','LineWidth',5)
% hold off
% saveX = getframe(ff).cdata;
% filename =[num2str(i),'_X.png'];
% imwrite(saveX ,filename);
xOld = sol1.value(x);
uOld = sol1.value(u);
x0new =xOld(:,timeshiftindx);
timeold=time1;
initialGuessX = [sol1.value(x(:,timeshiftindx:end)) repmat(sol1.value(x(:,end)),1,timeshiftindx-1)];
initialGuessU = [sol1.value(u(timeshiftindx:end)) repmat(sol1.value(u(end)),1,timeshiftindx-1)];
initialGuessX(6:7,:)=zeros(2,length(initialGuessX(1,:)));

end
disp('total time with warmstart was')
disp(num2str(toc))


%% New version MPC with pareto front plotting
% x_val  = [0 0 0 0 0 0 0]';
u_val = 0;
xf = [0,0]';

timehorizon = 20;
timestep = 0.5;
nHorizon = round(timehorizon/timestep);
nPoints = 20;
foh = true;
normalize = "fix";
% weight = [0.0514 0.9987];
weight = [0.5 0.5];

u_box = [-1, 1];


[ocp,x,u,d,x0,x0_p] = initializeOCP(timehorizon,timestep);
ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:(d.length()*timestep)]));
x_val = x0;
%%
% accum_costs = zeros(1,2);
costfun = [x(6,end), x(7,end)];

% ocp.set_value(x0_p, x_val(:,end))
% if foh
%     ocp.set_value(u(:,1), u_val(:,end))
% end
if isequal(normalize, "adapt")
    [p_params, up, np, norm_costfun, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize=normalize );
elseif isequal(normalize, "none")
    [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize=normalize );
elseif isequal(normalize, "fix")
    [p_params, ep_val, norm_costfun, ep, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize=normalize );
end
%%
sol = {};
J = {};
sol_init = [];
idx = uint8.empty(0,nHorizon);

for i = 1:nHorizon
    time = [(i-1)*timestep:timestep:timehorizon]
    ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep+((i-1)*timestep):timestep:(d.length()*timestep+((i-1)*timestep))]));
    if foh
        ocp.set_value(u(:,1), u_val(:,end))
    end
    if i > 1
         ocp.set_initial(x(:,2:end-1), sol{i-1}(idx(i-1)).value(x(:,3:end)));
         ocp.set_initial(u(:,1+foh:end-1), sol{i-1}(idx(i-1)).value(u(:,2+foh:end)));
        ocp.set_value(x(:,1), sol{i-1}(idx(i-1)).value(x(:,2)));
        ocp.set_value(u(:,1), sol{i-1}(idx(i-1)).value(u(:,2)));
    
    end

    if ismember(normalize, ["adapt", "none"])
        ep_fun_in = {ep_ocp, x0, x_val(:, end), w_ep, costfun};
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
    u_val(:, end+1) = sol{i}(idx(i)).value(u(:,2));

    ep_ocp.set_initial(x(:,2:end), sol_init.value(x(:,[3:end, end])))
    ep_ocp.set_initial(u(:,1+foh:end), sol_init.value(u(:,[2+foh:end, end])))
%     plot(x_val(1,:), x_val(2,:), sol{i}(idx(i)).value(x(1,2:end)), sol{i}(idx(i)).value(x(2,2:end)), "--", 'LineWidth', 2)
    time_val = ( [0:timestep:timestep*(length(x_val(1,:))-1)]);
time = [(i)*timestep:timestep: timestep*(length(sol{i}(idx(i)).value(x(1,2:end)))-1)+(i)*timestep];
plot(time_val, x_val(2,:), time, sol{i}(idx(i)).value(x(2,2:end)), "--", 'LineWidth', 2)
end

%%


function xdot = oszi_ode(x, u)
    xdot = [x(2); -x(1) + u; x(1:2)'*x(1:2); 0.5*u^2];
end

%%
function out = rescale_upnp(vec, up, np)

out = (vec - up)./(np - up);

end