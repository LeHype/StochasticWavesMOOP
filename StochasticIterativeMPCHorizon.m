% 1. Run the OPC with large time Horizon
%       - Find global weight for weighted sum. 
%       - 
global PathToParameters 
PathToParameters= 'src/PolySurge_inputs.mat';
load(PathToParameters);
filenameMOOP = ['MOOPStochastic_400seconds.mat'];
%%
timehorizon = 400;                          % shoud be self explanatory
timestep = 0.2;                             % shoud be self explanatory
nHorizon = round(timehorizon/timestep); 
                             % Number of Pareto Points
time = [0:timestep:(nHorizon)*timestep];    % Create array with discrete time steps
ocp_t = time;


[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(timehorizon,timestep);
monoW = monochromaticWave();
ocp.set_value(d,arrayfun(@(t) FBMStochasticWave(t),[0:timestep:((d.length()-1)*timestep)]));

costfun = ([x(6,end) x(7,end)]);
p_params = ocp.parameter(2,1);
ocp.minimize( costfun*p_params );
% [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method='ws', normalize='fix' );

tic
nPoints = 15;
w1= linspace(0.05, 1-0.05, 15);
w = [w1 ;1-w1]';
sol = [];
ResultsMOOP = cell(nPoints,1);
for i = 1:nPoints
    ocp.set_value(p_params,w(i,:));
    if ~isempty(sol)
        if i == 2
            solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
                'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
                'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
            ocp.solver('ipopt', struct(), solver_opt)
        end
%         ocp.set_initial(x(:,2:end),sol(i-1).value(x(:,2:end)))
        ocp.set_initial([ocp.x; ocp.lam_g], sol(i-1).value([ocp.x; ocp.lam_g]))
    end
    sol = [sol ocp.solve()];
    MOOPSolution = struct;
    MOOPSolution.weigths = w(i,:);
    MOOPSolution.x = ocp.value(x);
    MOOPSolution.u=ocp.value(u);
    MOOPSolution.ParetoPoint=[ocp.value(x(6,end)) ocp.value(x(7,end))];
    ResultsMOOP{i} = MOOPSolution;
end

% f1 = figure(1)
% 
% for i = 1:nPoints
%     scatter(sol(i).value(x(6,end)),sol(i).value(x(7,end)),250,'filled')
%     hold on
%     
% end
% l = legend();
% for i = 1:nPoints
% l.String{i} = num2str(i);
%     
% end

save(filenameMOOP,"ResultsMOOP");
toc
%% Chose Point here! Important choice. 
Point = 6;
% Validate with better integrator. 
d = sol(Point).value(d);


U = @(t) interp1(time,sol(Point).value(u),t,'previous');  %% Zero order hold
validationdgl = @(t,x) [Ac * x(1:5) - Bc * 1e6 * U(t) * gamma * x(2) + Bc * FBMStochasticWave(t)];
k = ode89(validationdgl,[0,timehorizon],x0(1:5));
% f2 = figure(2)
% plot(k.x,k.y(2,:))
% hold on 


% plot(time,sol(Point).value(x(2,:)))

RMSE = SignalDifference(k.x,k.y(2,:),time,sol(Point).value(x(2,:)))

%%
ocp_x = sol(Point).value(x);
ocp_u = sol(Point).value(u);
OCP = @(t) interp1(time,ocp_x(2,:),t);
OCP_U = @(t) interp1(time,ocp_u,t);
%% After confirming the time step is adequate start MPC here
Message = [];
for jj = 1:2
    
MPCtimehorizon = 60+5*jj;
wait_bar = waitbar(0,['Calculating MPC Horizon' num2str(MPCtimehorizon)]);
if (true)
filename = ['Stochastic,' datestr(now,'DD_HH_MM') '_MPC_Horizon_New_Step' num2str(MPCtimehorizon) 'Point_7.mat'];


% f3 = figure(3)
% subplot(3,1,1)
% plot(k.x,k.y(2,:),'r',LineWidth=5)
% hold on 
% subplot(3,1,2)
% plot(k.x,arrayfun(@(t)U(t),k.x),'r',LineWidth=5)
% hold on
% subplot(3,1,3)
% plot(k.x,arrayfun(@(t)U(t),k.x),'r',LineWidth=3)
% hold on


[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(MPCtimehorizon,timestep,x0=x0);
costfun = ([x(6,end) x(7,end)]);
p_params = ocp.parameter(2,1);

ocp.minimize( costfun*p_params );
ocp.set_value(p_params,w(Point,:));
starttime = [0:0.2:timehorizon-MPCtimehorizon];
ResultsMPC = cell(length(starttime),1);
AppliedSignal=[];
for i = 1:length(starttime)
    waitbar((i/length(starttime)),wait_bar,['Calculating MPC Horizon' num2str(MPCtimehorizon)]);
    nHorizon = round(MPCtimehorizon/timestep);     % Number of Pareto Points
    time = [starttime(i):timestep:(nHorizon)*timestep+starttime(i)];    % Create array with discrete time steps
    ocp.set_value(d,arrayfun(@(t) FBMStochasticWave(t),time));
    
    % if not the first MPC step
    if ~(i == 1)
        % reset solver for warmstarting
    ocp.set_value(x(:,1),solOld.value(x(:,2)));
    ocp.set_value(u(1),solOld.value(u(2)));

    ocp.set_initial(x(:,2:end-1),solOld.value(x(:,1:end-2)));
    ocp.set_initial(u(2:end-1),solOld.value(u(1:end-2)));
    end
    if i == 2
            solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
                'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
                'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
            ocp.solver('ipopt', struct(), solver_opt)
    end
    try
    solOld = ocp.solve();
    catch
    ocp.solver('ipopt')
 
    
%     Message = [Message ;'There was an error at increment' num2str(i) 'and MPCHorizon' num2str(MPCtimehorizon)]
    solOld = ocp.solve();

    solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
                'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
                'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
    ocp.solver('ipopt', struct(), solver_opt)

    end
    
    solMpc = struct;
    solMpc.x = solOld.value(x);
    solMpc.u = solOld.value(u);
    solMpc.time = solOld.value(time);
    AppliedSignal = [AppliedSignal [time(1); solOld.value(u(1))]];
    ResultsMPC{i} = solMpc;
%     subplot(3,1,1)
%     plot(time,solOld.value(x(2,:)))
%     hold on 
%     subplot(3,1,2)
%     plot(time,solOld.value(u))
%     hold on 
%     subplot(3,1,3)
%     plot(AppliedSignal(1,:),AppliedSignal(2,:))
%     hold on 
end

% 
% f4 = figure(4)
% for i = 1:length(starttime)
%   
% 
% Divergence = ResultsMPC{i}.u-arrayfun(@(t) OCP_U(t),ResultsMPC{i}.time);
% plot(ResultsMPC{i}.time,Divergence)
% hold on 

% f5 = figure(5)
% plot(AppliedSignal(2,:)-arrayfun(@(t) OCP_U(t),[0:timestep:AppliedSignal(1,end)]))

save(filename ,'ocp_t','OCP','OCP_U','ResultsMPC')
end
close(wait_bar)
end
%% TRY out different stuff to get turnpike visualized
% i = 10;
% f5 = figure(12)
% 
% plot(ResultsMPC{i}.time,(ResultsMPC{i}.u-arrayfun(@(t) OCP_U(t), ResultsMPC{i}.time)).^2)
% hold on 
% 









function [RMSE]  = SignalDifference(t1,sig1,t2,sig2)
sig2 = @(t) interp1(t2,sig2,t);
RMSE = 0;
for i = 1:length(t1)
    RMSE = RMSE + (sig1(i)-sig2(t1(i)))^2;



end
RMSE = sqrt(RMSE/length(t1));

if RMSE >= 0.02 
    warning('Please reconsider the time step! might be to large')
end
end