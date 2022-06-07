function [ocp,x,u,disturbance,x0] = initializeOCP(timehorizon,dt,params,args)
arguments
    timehorizon  (1,1) {mustBeNumeric}
    dt           (1,1) {mustBeNumeric}
    params       (1,:) {struct}
    args.solver  (1,:) {mustBeText} = 'ipopt'
    args.foh     (1,1) logical = true
    args.x0      (7,1) {mustBeNumeric} = zeros(7,1)
    args.ts      (1,1) {mustBeNumeric} = 0.0
end

%% Number of increments is timehorizon/dt (+1 if foh)
NumInc =round(timehorizon/(dt));

x0 = args.x0;
timeshift = args.ts;

%Construct the basic OCP from the ODE . only return the ode
%object. Will not apply disturbance and cost function. 





Path = zeros(200,8);



import casadi.*
load('/home/heib/Documents/HIWI/Flaßkamp/WaveHarvesting/DEA-Wave-Harvesting/PolySurge_inputs.mat', 'Ac', 'Bc', 'gamma');
slacks = MX;
wave_dgl = @(x,u,d) [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d;
                            cost_energy(x,u,d);
                            cost_damage(x,u)
                            ];


if all(x0 == 0) 
disp('I ran this')
    x0 = zeros(7,1);
    for i = 1:100
    if i <=100
    
        x0 = integrator_step_disturbed(x0,[0],1,wave_dgl,StochasticWave(i));
        Path(i,1:7) = transpose(full((evalf(x0))));
        Path(i,8) = StochasticWave(i);
    else
        x0= integrator_step_disturbed(x0,[0],1,wave_dgl,0);


Path(i,1:7) = transpose(full((evalf(x0))));
 
end
    end
x0=full(evalf(x0));
end
% To visualize execute:
if (false)
figure(1);
plot(Path(:,2));

hold on
yyaxis right
plot(Path(:,8))
end
% %%
% ode45(wave_dgl,10,)
%%

%% Basic implementation Just run once. 
timestep=dt;
x_box = [-Inf Inf; -Inf, Inf; -Inf Inf; -Inf, Inf; -Inf, Inf;-Inf, Inf;-Inf, Inf];
u_box = [0 33^2];

[ocp, x, u,varout{1:6}] = ode2ocp(wave_dgl, 7, 1, NumInc, dt, x0=x0, x_box=x_box, u_box=u_box, nd=1,foh=args.foh);
disturbance = varout{4};

ocp.set_value(disturbance,arrayfun(@(t) StochasticWave(t),[timestep+timeshift:timestep:(disturbance.length()*timestep)+timeshift]));

costs= ([-x(6,end) x(7,end)]-[params.up])./[params.np-params.up];
l   = ocp.variable( 1 );
pt  = ocp.parameter( 1, 2 );
dir = ocp.parameter( 1, 2 );
ocp.subject_to( pt + l*dir >= costs )      
ocp.minimize(-l)
ocp.set_value(pt,params.startingpoint);
ocp.set_value(dir,params.searchdirection);

%%
% DrawInput = @(i) plot(sol(i).value(u))
end
% figure
% for i = 1:16
%     subplot(4,4,i)
%     DrawInput(i)
% end

function x_end = integrator_step_disturbed(x0, u, dt, odefun, d)
% calculate one integration step with step size dt
import casadi.*

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
%%
% function xdot = wave_dgl(x,u,d)
% 
% 
% import casadi.*
% 
% load('/home/heib/Documents/HIWI/Flaßkamp/WaveHarvesting/DEA-Wave-Harvesting/PolySurge_inputs.mat', 'Ac', 'Bc', 'gamma');
% 
% slacks = MX;
% 
% % definition ode
% %  xdot = @(x, u, d, slacks) [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d;
% %                      cost_energy(x, u, d, slacks)];
%                    xdot =  @(x,u,d) Ac * x - Bc * 1e6 * u * gamma * x(2) + Bc * d;

% end


function cost = cost_energy(x, u, d)
    persistent Ch S R0
    
    if isempty(Ch) || isempty(S) || isempty(R0)
        load('/home/heib/Documents/HIWI/Flaßkamp/WaveHarvesting/DEA-Wave-Harvesting/PolySurge_inputs.mat', 'Ch', 'S', 'R0');
    end
    cost = (Ch*x(1).^2 + x(3:5)'*S*x(3:5) - d .* x(1))*1e-6 + u/R0;
end

function cost = cost_damage(x, u, ~, ~)   
%     cost = (max(u - 484/(cos(x(2)).^2), 0).^2)*1e-6;
    cost = (max(u - 484, 0).^2)*1e-6;
end

