global PathToParameters
PathToParameters= 'src/PolySurge_inputs.mat';
load(PathToParameters);             


%%
timehorizon     = 40;         % [1-inf]  How long
timestep        = 0.2;         % [0.05-1] MPC timestep, i.e. the discretisation of the ocp.
                               %          A step larger then 1 is not recommended.
SwingInTime     = 200;         % [100:~]  How long the system is left alone to swing in 
                               %          before the optimal control is applied. 
WaveForm        = 'Harmonic';  % ['Harmonic' 'Stochastic'] Choose the Wave Distrubance

saving          = true;       % If saving, the results will be saved to the "Results" folder.

filename = ['GroundTruth_Harmonic_SingleObjective.mat'];             
                               % If saving use this filename

nSteps          = round(timehorizon/timestep);       % Number of discrete timesteps
derivative_method = 'subgradient';


% Set the costfunction to add component to cost function simply set the 
% corresponting parameter to 0. For example params = [1 0 -1]' means 
% 1*first_component - 1*third_component
% i*u : params = [0 0 1 0 1 -1 -1 ]

params(1)   =   1;       % 1e-6*(Ch*x(1).^2)    
params(2)   =   1;       % 1e-6*x(3:5)'*S*x(3:5)
params(3)   =   1;       % u/R0
params(4)   =   -1;       % 1e-6*x(1)*d
params(5)   =   0;       % 0.5*C0*du
params(6)   =  -0;       % 0.5*gamma*x(2).^2*du
params(7)   =  -0;       % 2*gamma*x(1)*x(2)*u


 

%create OCP object and apply wave harvester DGL
[ocp,x,u,d,x0_p,du] = initializeOCPENERGY(timehorizon, timestep,ds=derivative_method,params=params);


% ocp.set_initial(x(:,2:end), warmstart.x(:,2:end))
% ocp.set_initial(u(2:end), warmstart.u(2:end))
ocp.solver('ipopt');
% Storage_Function = @(x,u)       0.5*Mh*x(1)^2 +0.5*Kh*x(2)^2+0.5*(C0-gamma*x(2)^2)*u + 0.5*x(3:5)'*Q*x(3:5); 

time            = linspace(0,timehorizon,d.length());% Create array with discrete time steps
WaveTime        = time+SwingInTime;                  % To create a smooth transition from the swing in the wave 
                                                     % continous after the swing in
                                                      
% ocp.subject_to(u(2:end)==0);
% Swing in the system for x seconds and set the initial value
x0 = SwingIn(SwingInTime, WaveForm, x0_p);
Storage_Function = @(x,u) 1e-6*(0.5*Mh*x(1)^2 +0.5*Kh*x(2)^2 + 0.5*x(3:5)'*Q*x(3:5))+ 0.5*(C0-gamma*x(2)^2)*u  -x(7); 

ocp.set_value(x0_p,x0);

% Apply the wave disturbance
switch WaveForm
    case 'Stochastic'
        ocp.set_value(d,arrayfun(@(t) FBMStochasticWave(t,Seed=Seed),WaveTime));
    case 'Harmonic'
        HarmonicWave = monochromaticWave();
        ocp.set_value(d,arrayfun(@(t) HarmonicWave(t),WaveTime));
end


costfun = (x(6,end)+ 0.2*x(7,end));

ocp.minimize(costfun);

ocp.solve()
sol = struct;
sol.x = ocp.value(x);
sol.u = ocp.value(u);
sol.time = time;
sol.d = ocp.value(d);
sol.du = ocp.value(du);
if (saving)
    if ~exist([pwd filesep 'Results'],'dir')
       mkdir('Results')
    end
    save(['Results' filesep filename],"sol")
end
%%
%Continue to TP_Zanon...
sol_old = sol;
%% Run validation ode89 simulation

 [x_ana,u_ana,d_ana,du_ana,t_ana]  = RunODE89Test(sol,params);
%%
%redefine parameters if you want
params(1)   =   -1;       % 1e-6*(Ch*x(1).^2)    
params(2)   =   -1;       % 1e-6*x(3:5)'*S*x(3:5)
params(3)   =   -1;       % u/R0
params(4)   =   +1;       % 1e-6*x(1)*d
params(5)   =   0;       % 0.5*C0*du
params(6)   =  -0;       % 0.5*gamma*x(2).^2*du
params(7)   =  -0;       % 2*gamma*x(1)*x(2)*u
%%

SF = [];

cost_from_ocp = [sol.x(6,1) sol.x(6,2:end)-sol.x(6,1:end-1)];  % the cost that was used in ocp
du = [ 0 ((sol.u(3:end)-sol.u(2:end-1))/timestep +  (sol.u(2:end-1)-sol.u(1:end-2))/timestep)/2  ((sol.u(end)-sol.u(end-1))/timestep +  (sol.u(end-1)-sol.u(end-2))/timestep)/2 ];

cost_recalculated= [];                              % maybe if we want to calcualte other costs after the fact
                                % If you want to change the SR 

for i = 1:length(sol.x)
    SF = [SF Storage_Function(sol.x(:,i),sol.u(i))];
    cost_recalculated = [cost_recalculated cost_energy(sol.x(:,i),sol.u(i),sol.d(i),du(i),params)];
end

SF_ana = [];


cost_from_ocp_ana = [x_ana(6,1) x_ana(6,2:end)-x_ana(6,1:end-1)];  % the cost that was used in ocp

cost_recalculated_ana= [];                              % maybe if we want to calcualte other costs after the fact
                                    % If you want to change the SR 
for i = 1:length(x_ana)
    SF_ana = [SF_ana Storage_Function(x_ana(:,i),u_ana(i))];
    cost_recalculated_ana = [cost_recalculated_ana cost_energy(x_ana(:,i),u_ana(i),d_ana(i),du_ana(i),params)];
end


xlimit = [100 120];
figure(2)

plot(t_ana(2:end),SF_ana(2:end)-SF_ana(1:end-1))
hold on
% plot(t_ana,cost_recalculated_ana)
plot(t_ana,cost_from_ocp_ana)
title('Storage function and Supply rate, Analytical')
xlim (xlimit)
legend('[Ana] Difference of storage function', '[Ana] Supply rate')
EGFixFigure

figure (3)
plot(time(2:end),SF(2:end)-SF(1:end-1))
hold on
plot(time,cost_from_ocp)

xlim (xlimit)



title('Storage function and Supply rate from OCP')
legend('[OCP] Difference of storage function', '[OCP] supply rate ocp')
ylabel ('??')

EGFixFigure
%%
figure(5)
plot(t_ana,u_ana)
hold on 
plot(time,sol.u)
xlabel('time (s)')
ylabel('voltage squared (v^2)')
xlim (xlimit)
legend('voltage analytical' , 'voltage ocp')
EGFixFigure

figure(6)
plot(t_ana,rad2deg(x_ana(1,:)))
hold on 
plot(time,rad2deg(sol.x(1,:)))
xlabel('time (s)')
ylabel('angle of flap (°)')
xlim (xlimit)
legend('analytical angle ' , 'angle from ocp')
EGFixFigure
%%
plot(t_ana(2:end),SF_ana(2:end)-SF_ana(1:end-1)-cost_from_ocp_ana(2:end) )
yline(0)

title('reconst')