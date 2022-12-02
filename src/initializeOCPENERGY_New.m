function [ocp,x,u,disturbance,x0_p,du] = initializeOCPENERGY(timehorizon,dt,args)
arguments
    timehorizon     (1,1) {mustBeNumeric}
    dt              (1,1) {mustBeNumeric}
    args.solver     (1,:) {mustBeText} = 'ipopt'
    args.foh        (1,1) logical = true
    args.params     (:,:) {mustBeNumeric} = [0 0 1 0 1 -1 -1 ]
    args.ds      (1,:) {mustBeText} = 'central'     %derivative scheme
end


global PathToParameters %% Just so one does not have to change it in multiple scripts
%% Number of increments is timehorizon/dt (+1 if foh)
NumInc =round(timehorizon/(dt));



%Construct the basic OCP from the ODE . only return the ode
%object. Will not apply disturbance and cost function. 
% System State is denoted by x= [theta_dot theta z]' 
load(PathToParameters ,'Ac', 'Bc', 'gamma');

wave_dgl = @(x,u,d,du) [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d;
                        cost_energy(x,u,d,du,args.params);
                        d*x(1)*1e-6
                        ];
nx = 7;




%% Basic implementation  
u_box = [0 (33^2)];

[ocp, x, u,varout{1:8}] = ode2ocp_new(wave_dgl, nx, 1, NumInc, dt, x0='param', u_box=u_box, nd=1, foh=args.foh,ds=args.ds);
du = varout{8};
ocp.set_value(u(1),0);       %% first value for u has to be zero. 

% ocp.set_value(du(1),0);
disturbance = varout{4};     
x0_p = varout{2};            %% Need to include the object of x0 so i can change it outside of the function


end






