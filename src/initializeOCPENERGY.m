function [ocp,x,u,disturbance,x0,x0_p] = initializeOCP(timehorizon,dt,args)
arguments
    timehorizon  (1,1) {mustBeNumeric}
    dt           (1,1) {mustBeNumeric}
    args.solver  (1,:) {mustBeText} = 'ipopt'
    args.foh     (1,1) logical = true
    args.x0      (12,1) {mustBeNumeric} = zeros(12,1)
end


global PathToParameters %% Just so one does not have to change it in multiple scripts
%% Number of increments is timehorizon/dt (+1 if foh)
NumInc =round(timehorizon/(dt));

x0 = args.x0;

%Construct the basic OCP from the ODE . only return the ode
%object. Will not apply disturbance and cost function. 

import casadi.*
load(PathToParameters ,'Ac', 'Bc', 'gamma','R0');
% slacks = MX;
wave_dgl = @(x,u,d) [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d;
                            cost_energy(x,u,d);
                            cost_damage(x,u);
                            EE1(x,u,d);
                            EE2(x,u,d);
                            EE3(x,u,d);
                            EE4(x,u,d);
                            gamma*x(2)*x(1)*u-u/(R0)
                            ];

%If no x0 is given run a simulation for x seconds so swing
% the system. The behavior of the system may vary for different
% initial states but only for the first couple of seconds

if all(x0 == 0) 

    x0 = zeros(12,1);
    for i = 1:100  %%100 seconds swing in 

    
        x0 = integrator_step_disturbed(x0,[0],1,wave_dgl,NewStochasticWave(i));





 

    end
x0=full(evalf(x0));
end
x0(6)=0;
x0(7)=0;
%%make sure the energy is reset to zero
% To visualize execute:
%%

%% Basic implementation  
x_box = [-Inf Inf; -Inf, Inf; -Inf Inf; -Inf, Inf; -Inf, Inf;-Inf, Inf;-Inf, Inf;-Inf, Inf;-Inf, Inf;-Inf, Inf;-Inf, Inf;-Inf, Inf];
u_box = [0 33^2];

[ocp, x, u,varout{1:6}] = ode2ocp(wave_dgl, 12, 1, NumInc, dt,x0='no', x_box=x_box, u_box=u_box, nd=1,foh=args.foh);
ocp.set_value(u(1),0);       %% first value for u has to be zero. 
disturbance = varout{4};     
x0_p = varout{2};            %% Need to include the object of x0 so i can change it outside of the function

ocp.set_value(x0_p,x0);



end


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
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters ,'Ch', 'S', 'R0');
    end
    cost = (Ch*x(1).^2 + x(3:5)'*S*x(3:5) - d .* x(1))*1e-6 + u/R0;
end
function E1 = EE1(x, u, d)
    persistent Ch S R0
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters, 'Ch', 'S', 'R0');
    end
    E1 = 1e-6*(Ch*x(1).^2);
end
function E2 = EE2(x, u, d)
    persistent Ch S R0
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters, 'Ch', 'S', 'R0');
    end
     E2 = (x(3:5)'*S*x(3:5))*1e-6 ;
end
function E3 = EE3(x, u, d)
    persistent Ch S R0
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters, 'Ch', 'S', 'R0');
    end
      E3 = - (d .* x(1))*1e-6;
end
function E4 = EE4(x, u, d)
    persistent Ch S R0
    global PathToParameters
    if isempty(Ch) || isempty(S) || isempty(R0)
        load(PathToParameters, 'Ch', 'S', 'R0');
    end
     E4 = u/R0;
end

function cost = cost_damage(x, u, ~, ~)   
%     cost = (max(u - 484/(cos(x(2)).^2), 0).^2)*1e-6;
    cost = (max(u - 484, 0).^2)*1e-6;
end
