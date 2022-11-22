global PathToParameters
PathToParameters= 'src/PolySurge_inputs.mat';
load(PathToParameters);             

load("tp_explicit_1_0.mat");

%%
timehorizon     = 200;         % [1-inf]  How long
timestep        = 0.2;         % [0.05-1] MPC timestep, i.e. the discretisation of the ocp.
                               %          A step larger then 1 is not recommended.
SwingInTime     = 200;         % [100:~]  How long the system is left alone to swing in 
                               %          before the optimal control is applied. 
WaveForm        = 'Harmonic';  % ['Harmonic' 'Stochastic'] Choose the Wave Distrubance

saving          = false;       % If saving, the results will be saved to the "Results" folder.

filenameMOOP = ['MOOPStochastic_400seconds.mat'];             
                               % If saving use this filename

nSteps          = round(timehorizon/timestep);       % Number of discrete timesteps

%create OCP object and apply wave harvester DGL
[ocp,x,u,d,x0_p] = initializeOCPENERGY(timehorizon, timestep, get_energy=true);
ocp.solver('ipopt');
Storage_Function = @(x,u) 1e-6*(0.5*Mh*x(1)^2 +0.5*Kh*x(2)^2 + 0.5*x(3:5)'*Q*x(3:5)) + 0.5*(C0-gamma*x(2)^2)*u; 

time            = linspace(0,timehorizon,d.length());% Create array with discrete time steps
WaveTime        = time+SwingInTime;                  % To create a smooth transition from the swing in the wave 
                                                     % continous after the swing in
                                                      
% ocp.subject_to(u(2:end)==0);
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
SR_sym = (x(6,2:end) - x(6,1:end-1)) - repmat(tp.ell, 1, (size(x,2)-1)/length(tp.ell));
costfun = sum(SR_sym);

ocp.minimize(-costfun);

ocp.solve()
SR = ocp.value(SR_sym);
sol = struct;
sol.x = ocp.value(x);
sol.u = ocp.value(u);
sol.time = time;
sol.d = ocp.value(d);
%%
SF = [];
cost_implicit = [];
cost_explicit = [];
x_rep = repmat(tp.x, 1, (size(x,2)-1)/length(tp.ell));
u_rep = repmat(tp.u, 1, (length(u)-1)/length(tp.ell));
for i = 1:length(sol.x)-1
    SF = [SF Storage_Function(sol.x(1:5,i) - x_rep(:,i), sol.u(i) - u_rep(:,i))];
end
%%
figure(1)

plot(time,sol.x(6,:))
hold on 
plot(time,-sol.x(12,:))
xlabel('time [s]')
ylabel('Energy Harvested')
EGFixFigure;
l = legend('implicit energy calculation $\int_{0}^{t} C_h \dot{\theta}^2+x^T S x+\frac{u}{R_0}-d(t) \dot{\theta}$' , 'explicit energy calculation $(\int_{0}^{t} \gamma*\theta*\dot{\theta}u-u/R_0$');
set(l,'Interpreter','Latex');


figure(2)
plot(time,-cost_implicit)
hold on 
plot(time,cost_explicit)
xlabel('time [s]')
ylabel('Energy Harvested')
EGFixFigure;
l = legend('implicit energy calculation $-C_h \dot{\theta}^2+x^T S x+\frac{u}{R_0}-d(t) \dot{\theta}$' , 'explicit energy calculation $(\gamma*\theta*\dot{\theta}u-u/R_0$');

set(l,'Interpreter','Latex');


ff = figure(3);
title("Energy decomposition");

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

