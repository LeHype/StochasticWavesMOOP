% 1. Run the OPC with large time Horizon
%       - Find global weight for weighted sum. 
%       - 
global PathToParameters 
PathToParameters= '/home/heib/Documents/HIWI/Flaßkamp/WaveHarvesting/DEA-Wave-Harvesting/PolySurge_inputs.mat';
load(PathToParameters);
%%
timehorizon = 100;                           % shoud be self explanatory
timestep = 0.2;                             % shoud be self explanatory
nHorizon = round(timehorizon/timestep); 
                             % Number of Pareto Points
time = [0:timestep:(nHorizon)*timestep];    % Create array with discrete time steps
MPCtimehorizon = 10;
filename = ['Harmonic_Wave,' datestr(now,'DD_HH_MM') '.mat'];

[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(timehorizon,timestep);
monoW = monochromaticWave();
ocp.set_value(d,arrayfun(@(t) monoW(t),[0:timestep:((d.length()-1)*timestep)]));

costfun = ([x(6,end) x(7,end)]);
p_params = ocp.parameter(2,1);
ocp.minimize( costfun*p_params );
% [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method='ws', normalize='fix' );
%%
tic
nPoints = 10;
w1= [0.05:1/(nPoints):1-0.05];
w = [w1 ;1-w1]';
sol = [];
for i = 1:nPoints
    ocp.set_value(p_params,w(i,:));
    if ~isempty(sol)
         ocp.set_initial(x(:,2:end),sol(i-1).value(x(:,2:end)))
    end
    sol = [sol ocp.solve()];
end

f1 = figure(1)

for i = 1:nPoints
    scatter(sol(i).value(x(6,end)),sol(i).value(x(7,end)),250,'filled')
    hold on
    
end
l = legend();
for i = 1:nPoints
l.String{i} = num2str(i);
    
end


toc
%% Chose Point here! Important choice. 
Point = 5;
% Validate with better integrator. 
d = sol(Point).value(d);


U = @(t) interp1(time,sol(Point).value(u),t,'previous');  %% Zero order hold
validationdgl = @(t,x) [Ac * x(1:5) - Bc * 1e6 * U(t) * gamma * x(2) + Bc * monoW(t)];
k = ode89(validationdgl,[0,100],x0(1:5));
f2 = figure(2)
plot(k.x,k.y(2,:))
hold on 

plot(time,sol(Point).value(x(2,:)))

RMSE = SignalDifference(k.x,k.y(2,:),time,sol(Point).value(x(2,:)))

OPC = @(t) interp1(time,sol(Point).value(x(2,:)),t);
OPC_U = @(t) interp1(time,sol(Point).value(u),t);
%% After confirming the time step is adequate start MPC here

f3 = figure(3)
subplot(3,1,1)
plot(k.x,k.y(2,:),'r',LineWidth=5)
hold on 
subplot(3,1,2)
plot(k.x,arrayfun(@(t)U(t),k.x),'r',LineWidth=5)
hold on
subplot(3,1,3)
plot(k.x,arrayfun(@(t)U(t),k.x),'r',LineWidth=3)
hold on


[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(MPCtimehorizon,timestep,x0=x0);
costfun = ([x(6,end) x(7,end)]);
p_params = ocp.parameter(2,1);

ocp.minimize( costfun*p_params );
ocp.set_value(p_params,w(Point,:));
starttime = [0:1:timehorizon-MPCtimehorizon];
ResultsMPC = cell(length(starttime),1);
AppliedSignal=[];
for i = 1:length(starttime)
    nHorizon = round(MPCtimehorizon/timestep);     % Number of Pareto Points
    time = [starttime(i):timestep:(nHorizon)*timestep+starttime(i)];    % Create array with discrete time steps
    ocp.set_value(d,arrayfun(@(t) monoW(t),time));
    
    % if not the first MPC step
    if ~(i == 1)
    ocp.set_value(x(:,1),solOld.value(x(:,6)));
    ocp.set_value(u(1),solOld.value(u(6)));

    ocp.set_initial(x(:,2:end-1),solOld.value(x(:,1:end-2)));
    ocp.set_initial(u(2:end-1),solOld.value(u(1:end-2)));
    end
    
    solOld = ocp.solve();
    solMpc = struct;
    solMpc.x = solOld.value(x);
    solMpc.u = solOld.value(u);
    solMpc.time = solOld.value(time);
    AppliedSignal = [AppliedSignal [time(1:5); solOld.value(u(1:5))]];
    ResultsMPC{i} = solMpc;
    subplot(3,1,1)
    plot(time,solOld.value(x(2,:)))
    hold on 
    subplot(3,1,2)
    plot(time,solOld.value(u))
    hold on 
    subplot(3,1,3)
    plot(AppliedSignal(1,:),AppliedSignal(2,:))
    hold on 
end


%% Plot the difference to the MPC over time
f4 = figure(4)
for i = 1:length(starttime)
  
Divergence = ResultsMPC{i}.u-arrayfun(@(t) OPC_U(t),ResultsMPC{i}.time);
plot(ResultsMPC{i}.time,Divergence)
hold on 
end
f5 = figure(5)
plot(AppliedSignal(2,:)-arrayfun(@(t) OPC_U(t),[0:timestep:AppliedSignal(1,end)]))




function [RMSE]  = SignalDifference(t1,sig1,t2,sig2)
sig2 = @(t) interp1(t2,sig2,t);
RMSE = 0;
for i = 1:length(t1)
    RMSE = RMSE + (sig1(i)-sig2(t1(i)))^2;
save(filename ,'OPC','OPC_U','ResultsMPC','f1','f2','f3','f4','f5')

end
RMSE = sqrt(RMSE/length(t1));

if RMSE >= 0.02 
    warning('Please reconsider the time step! might be to large')
end
end