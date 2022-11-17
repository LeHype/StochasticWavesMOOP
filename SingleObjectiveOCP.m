global PathToParameters
PathToParameters= 'src/PolySurge_inputs.mat';
load(PathToParameters);



%%
timehorizon     = 300;         % [1-inf]  How long
timestep        = 0.5;         % [0.05-1] MPC timestep, i.e. the discretisation of the ocp.
                               %          A step larger then 1 is not recommended.
SwingInTime     = 200;         % [100:~]  How long the system is left alone to swing in 
                               %          before the optimal control is applied. 
WaveForm        = 'Harmonic';  % ['Harmonic' 'Stochastic'] Choose the Wave Distrubance

saving          = false;       % If saving, the results will be saved to the "Results" folder.

filenameMOOP = ['MOOPStochastic_400seconds.mat'];             
                               % If saving use this filename

nSteps          = round(timehorizon/timestep);       % Number of discrete timesteps

%create OCP object and apply wave harvester DGL
[ocp,x,u,d,x0_p] = initializeOCPENERGY(timehorizon,timestep);
ocp.solver('ipopt');
Storage_Function = @(x,u) 0.5*Mh*x(1)^2 +0.5*Kh*x(2)^2+0.5*(C0-gamma*x(2)^2)*u +0.5*x(3:5)'*Q*x(3:5); 

time            = linspace(0,timehorizon,d.length());% Create array with discrete time steps
WaveTime        = time+SwingInTime;                  % To create a smooth transition from the swing in the wave 
                                                     % continous after the swing in
                                                      

% Swing in the system for x seconds and set the initial value
x0 = SwingIn(SwingInTime, WaveForm, x0_p);

ocp.set_value(x0_p,x0);

% Apply the wave disturbance
switch WaveForm
    case 'Stochastic'
        ocp.set_value(d,arrayfun(@(t) FBMStochasticWave(t,Seed=Seed),WaveTime));
    case 'Harmonic'
        HarmonicWave = monochromaticWave();
        ocp.set_value(d,arrayfun(@(t) HarmonicWave(t),WaveTime));
end
% Set the costfunction
costfun = (x(6,end));

ocp.minimize(costfun);

ocp.solve()
sol = struct;
sol.x = ocp.value(x);
sol.u = ocp.value(u);
sol.time = time;
sol.d = ocp.value(d);
%%
SF = [];
cost_implicit = [];
cost_explicit = [];
for i = 1:length(sol.x)
    SF = [SF Storage_Function(sol.x(:,i),sol.u(i))];
    cost_implicit = [cost_implicit cost_energy(sol.x(:,i),sol.u(i),sol.d(i))];
    cost_explicit = [cost_explicit cost_energy_explicit(sol.x(:,i),sol.u(i),sol.d(i))];

end
%%
figure(1)

plot(time,sol.x(6,:))
hold on 
plot(time,-sol.x(12,:))
xlabel('time [s]')
ylabel('Energy Harvested')
EGFixFigure;
l = legend('inplicit energy calculation $\int_{0}^{t} C_h \dot{\theta}^2+x^T S x+\frac{u}{R_0}-d(t) \dot{\theta}$' , 'explicit energy calculation $(\int_{0}^{t} \gamma*\theta*\dot{\theta}u-u/R_0$');
set(l,'Interpreter','Latex');


figure(2)
plot(time,-cost_implicit)
hold on 
plot(time,cost_explicit)
xlabel('time [s]')
ylabel('Energy Harvested')
EGFixFigure;
l = legend('inplicit energy calculation $-C_h \dot{\theta}^2+x^T S x+\frac{u}{R_0}-d(t) \dot{\theta}$' , 'explicit energy calculation $(\gamma*\theta*\dot{\theta}u-u/R_0$')

set(l,'Interpreter','Latex');


ff = figure(3)
title("Energy decomposition")

plot(time,sol.x(6,:));
hold on 
plot(time,sol.x(8,:));
plot(time,sol.x(9,:));
plot(time,sol.x(10,:));
plot(time,sol.x(11,:));

% plot(time,Validation,'--');
plot(time,-sol.x(12,:));
xlabel(['Time'])
legend(['Total Energy'],['$C_h \dot\Theta^2$'],['$\frac{1}{2} z^T S_r z$'],['$-d \dot\Theta$'],['$\frac{u}{R_0}$'],['Direct Energy Calculation'],'Interpreter','latex','Fontsize',22)

% p_params = ocp.parameter(2,1);
% ocp.minimize( costfun*p_params );
% % [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method='ws', normalize='fix' );
% 
% tic
% nPoints = 15;
% w1= linspace(0.05, 1-0.05, 15);
% w = [w1 ;1-w1]';
% sol = [];
% ResultsMOOP = cell(nPoints,1);
% for i = 1:nPoints
%     ocp.set_value(p_params,w(i,:));
%     if ~isempty(sol)
%         if i == 2
%             solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
%                 'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
%                 'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
%             ocp.solver('ipopt', struct(), solver_opt)
%         end
%        ocp.set_initial([ocp.x; ocp.lam_g], sol(i-1).value([ocp.x; ocp.lam_g]))
%     end
%     sol = [sol ocp.solve()];
%     MOOPSolution = struct;
%     MOOPSolution.weigths = w(i,:);
%     MOOPSolution.x = ocp.value(x);
%     MOOPSolution.u=ocp.value(u);
%     MOOPSolution.ParetoPoint=[ocp.value(x(6,end)) ocp.value(x(7,end))];
%     ResultsMOOP{i} = MOOPSolution;
% end
function cost = cost_energy(x, u, d)
    persistent Ch S R0
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters ,'Ch', 'S', 'R0');
    end
    cost = (Ch*x(1).^2 + x(3:5)'*S*x(3:5) - d .* x(1))*1e-6 + u/R0;
end
function cost = cost_energy_explicit(x, u, d)
    persistent Ch S R0 gamma
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters ,'Ch', 'S', 'R0','gamma');
    end
    cost = gamma*x(2)*x(1)*u-u/(R0);
end

