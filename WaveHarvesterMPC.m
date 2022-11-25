%% Runs MPC algorithm for wave harvester
% This script requires the installation of the CasADi toolbox, available at https://web.casadi.org/

% Define MPC parameters. These can and should be chosen by the user.
                               
                               % Suggested
                               % values

MPCtimehorizon  = 60;          % [1-70]   MCP Timehorizon. After ~70 seconds the solution does not improve.
simulation_time = 3000;        % [1-inf]  How long
timestep        = 0.5;         % [0.05-1] MPC timestep, i.e. the discretisation of the ocp.
                               %          A step larger then 1 is not recommended.
Seed            = 4  ;         % [1-10]   Seed of the Wave distrurbance. Seeds [1-10] have been provided.
Damagereduction = 0.4;         % [0-1]    This implementation uses multi-criteria Optimization.
                               %          Therefore we need to specify a weight for the damage contribution. 
                               %          0 would correspond to single-objective optimization only
                               %          maximizing the harvested energy
SwingInTime     = 200;         % [100:~]  How long the system is left alone to swing in 
                               %          before the optimal control is applied. 
WaveForm        = 'Stochastic';% ['Harmonic' 'Stochastic'] Choose the Wave Distrubance
excitation_factor = 1;         % [0.1-2]  Factor for s aling the wave excitation. Warning: The model's accuracy
                               %          is reduced for theta greater than 30°.
                                      

plotting        = true;        % If plotting the solution of the MPC is shown as it progresses.
                               % but slows down the algorithm!
saving          = true;        % If saving, the results will be saved to the "Results" folder.
filename = ['MPC' WaveForm '_' num2str(simulation_time) 'seconds_Damagereduction_' num2str(Damagereduction) '_Seed_' num2str(Seed) '.mat'];             
                               % If saving use this filename             

% Here, the algorithm starts. User interference is not recomended


weight = [1-Damagereduction Damagereduction];       % Weight for Weighted Sum MOOP

[ocp,x,u,d,x0_p] = initializeOCPENERGY(MPCtimehorizon,timestep);                                             

% Swing in the system for x seconds and set the initial value
x0 = SwingIn(SwingInTime,WaveForm,Seed=Seed);
ocp.set_value(x0_p,x0);

% Set the costfunction and weights
costfun = ([x(6,end), x(7,end)]);
p_params = ocp.parameter(2, 1);
ocp.minimize( costfun*p_params );
ocp.set_value(p_params, weight);
% Apply the wave disturbance
switch WaveForm
    case 'Stochastic'
        Wave = @(t) FBMStochasticWave(t,Seed=Seed);
        
    case 'Harmonic'
        Wave = monochromaticWave();
end

% the MPC starts here. The MPC solutions are tracked in the
% ResultsMPC array. The applied signal is stored in the AppliedSignal array
% AppliedSignal (1,:) : Global time
% AppliedSignal (2,:) : Applied Voltage
% AppliedSignal (3,:) : Angle of Flap

starttime = 0+SwingInTime:timestep:simulation_time+SwingInTime;
ResultsMPC = cell(length(starttime),1);
AppliedSignal=[];
if (plotting)
    close all
    f = openfig(['src' filesep 'MPCPlotting.fig']);
    ax1 = f.Children(3);
    ax2 = f.Children(2);
    ax3 = f.Children(1);
end

for i = 1:length(starttime)
    nHorizon = round(MPCtimehorizon/timestep);                          % Number of Pareto Points
    time = [starttime(i):timestep:(nHorizon)*timestep+starttime(i)];    % Create array with discrete time steps
    ocp.set_value(d,arrayfun(@(t) Wave(t),time));

    % if not the first MPC step set u_0 and x_0 to the solution of
    % previous MPC.
    if ~(i == 1)
        % reset solver for warmstarting

        ocp.set_value(x(:,1), solOld.value(x(:,2)));
        ocp.set_value(u(1),   solOld.value(u(2)));
        
        ocp.set_initial(ocp.x,     solOld.value(ocp.x));
        ocp.set_initial(ocp.lam_g, solOld.value(ocp.lam_g));
    end

    if i == 2
        % Set the solver to warmstart
        %  (faster because we have good initial guess from precious MPC)
        solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
            'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
            'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8, 'print_level', 0);
        ocp.solver('ipopt', struct(), solver_opt)
    end
    try
        solOld = ocp.solve();
    catch
        % The chosen warmstart parameters do not always work. If that is the case
        % disable warmstart for one iteration step.
        ocp.solver('ipopt', struct, struct('print_level', 0))
        solOld = ocp.solve();
        solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
            'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
            'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8, 'print_level', 0);
        ocp.solver('ipopt', struct(), solver_opt)

    end

    %Save results in struct
    solMpc = struct;
    solMpc.x = solOld.value(x);
    solMpc.u = solOld.value(u);
    solMpc.time = solOld.value(time);
    AppliedSignal = [AppliedSignal [time(1); solOld.value(u(1)) ; solOld.value(x(2,1))]];
    ResultsMPC{i} = solMpc;


    if (plotting)
        axes(ax1);
        hold on
        plot(AppliedSignal(1,:),rad2deg(AppliedSignal(3,:)),'b',LineWidth=4)
        plot(solOld.value(time),rad2deg(solOld.value(x(2,:))),'r--',LineWidth=4)
        
        if ~(i==1)
            delete(ax1.Children(end))
            delete(ax1.Children(end))
        end

        axes(ax2);
        hold on
        plot(AppliedSignal(1,:),AppliedSignal(2,:),'b',LineWidth=4)
        plot(solOld.value(time),solOld.value(u),'r--',LineWidth=4)
        
        if ~(i==1)
            delete(ax2.Children(end))
            delete(ax2.Children(end))
        end
        
        axes(ax3);
        hold on
        plot(AppliedSignal(1,:),arrayfun(@(t) Wave(t),AppliedSignal(1,:)),'b',LineWidth=4)
        plot(solOld.value(time),arrayfun(@(t) Wave(t),solOld.value(time)),'r--',LineWidth=4)
        
        if ~(i==1)
            delete(ax3.Children(end))
            delete(ax3.Children(end))
        end
    end

end
if (saving)
    if ~exist([pwd filesep 'Results'],'dir')
       mkdir('Results')
    end
    save(['Results' filesep filename],"ResultsMPC",'time')
end

%% Damage cost function
function cost = cost_damage(x, u, ~, ~)
%     cost = (max(u - 484/(cos(x(2)).^2), 0).^2)*1e-6;
cost = (max(u - 484, 0).^2)*1e-6;
end

%% Energy cost function
function cost = cost_energy(x, u, d)
persistent Ch S R0

if isempty(Ch) || isempty(S) || isempty(R0)
    load(['src' filesep 'PolySurge_inputs.mat'] ,'Ch', 'S', 'R0');
end
cost = (Ch*x(1).^2 + x(3:5)'*S*x(3:5) - d .* x(1))*1e-6 + u/R0;
end


%% runge-kutta 4  integrator
function x_end = integrator_step_disturbed(x0, u, dt, odefun, d)
% calculate one integration step with step size dt

x0_rk = x0;
k = casadi.MX( size( x0, 1 ), 4 );

if size(u,2) == 1
    k(:,1) = odefun(x0_rk(:,end)                  , u, d);
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u, d);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u, d);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u, d);
else
    k(:,1) = odefun(x0_rk(:,end)                  , u(:,1),       d(:,1));
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u(:,2),       d(:,2));
end

x_end  = x0_rk + dt / 6 * k * [1 2 2 1]';

end

