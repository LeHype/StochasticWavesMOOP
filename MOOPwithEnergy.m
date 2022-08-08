%% This is an example usage for the multiobjective OCP
% In this formulation the individual energy terms are part
% of the calculation. The state vector looks as follows:
%x1-x5 state variables in given order (like q in paper)
%x6 and x7 are the total energy (integral) and damage
%x8, x9, x10, and x11 are the energy contributions in the order they appear
%%
% Adjustment Needed! Change path here. 
global PathToParameters 
PathToParameters= '/home/heib/Documents/HIWI/Flaßkamp/WaveHarvesting/DEA-Wave-Harvesting/PolySurge_inputs.mat';
load(PathToParameters);





%% 
% The initial state of the system. It is not usefull to 
% start the system at [ 0 0 .. 0 0]. It's better to have it already
% somewhat swung up. 
% To avoid the simulation for the swing up uncomment the next line.

% x0 = [      0.0961 -0.2059 214.7006     293.8507    261.3057 0 0 0 0 0 0 0]';

timehorizon = 40;                           % shoud be self explanatory
timestep = 0.7;                             % shoud be self explanatory
nHorizon = round(timehorizon/timestep); 
nPoints = 6;                                % Number of Pareto Points
time = [0:timestep:(nHorizon)*timestep];    % Create array with discrete time steps

% the creation of the ocp object is encapsulated in this function.
% x,u, ect need to be returned to access these objects later
[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(timehorizon,timestep);


%Set the disturbance to the wave function
ocp.set_value(d,arrayfun(@(t)NewStochasticWave(t),[0:timestep:((d.length()-1)*timestep)]));

%Define Storage Function (Not relevant for OCP)
Storage_Function =  @(x,u) 0.5*Mh*x(1)^2 +0.5*Kh*x(2)^2+0.5*(C0-gamma*x(2)^2)*u +0.5*x(3:5)'*Q*x(3:5); 

% Cost is energy x6 and damage x7
costfun = ([x(6,end) x(7,end)]);
% Scalarize MOOCP just so for example a weight of 0.5 
% always represents the point in the middle of the 
% Paretofront
[p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize='fix' );

%% Calculate a single Point on the Pareto front 
POINT = 1;   
ep = ep_ocp;
weights = linspace(1,0,nPoints);
weights(1) = weights(1)-1E-3;
weights(end) = weights(end) + 1E-3;
weights = [weights; 1-weights];
[spt, sdir] = p_params{:};              % get starting point and search direction parameters

ndir = null(ep(2:end,:)-ep(1,:));
if all(ndir >= 0)
    ndir = -ndir;
end
ocp.set_value(sdir, ndir);

bp_pts = weights(1,:)'.*ep(1,:) + weights(2,:)'.*ep(2,:);


    ocp.set_value( spt, bp_pts(POINT, :) )
    sol = ocp.solve();



Validation = [sol.value(x(8,:))+sol.value(x(9,:))+sol.value(x(10,:))+sol.value(x(11,:))];
%% Just plotting here. The Energy components are only for visulization and not for the OCP

f = figure(3)

subplot(1,3,1)
plot(time, sol.value(x(2,:)))
xlim([0 time(end)])
f.CurrentAxes.FontSize = 25;
ylabel('Angle $\theta$','Interpreter','latex')
xlabel('Time')
subplot(1,3,2)
plot(time, sol.value(u))
xlim([0 time(end)])
xlabel('Time')
ylabel('Input U, $V^2$','Interpreter','latex')
f.CurrentAxes.FontSize = 25;
subplot(1,3,3)
plot(time, sol.value(d))
xlim([0 time(end)])
ylabel('Wave disturbance')
xlabel('Time')
f.CurrentAxes.FontSize = 25;
set(f,'color','w');
%%
ff = figure(2)
title("Energy decomposition")
set(0, 'DefaultLineLineWidth', 4);
plot(time,sol.value(x(6,:)));
hold on 
plot(time,sol.value(x(8,:)));
plot(time,sol.value(x(9,:)));
plot(time,sol.value(x(10,:)));
plot(time,sol.value(x(11,:)));

% plot(time,Validation,'--');
plot(time,-sol.value(x(12,:)));
xlabel(['Time'])
legend(['Total Energy'],['$C_h \dot\Theta^2$'],['$\frac{1}{2} z^T S_r z$'],['$-d \dot\Theta$'],['$\frac{u}{R_0}$'],['Direct Energy Calculation'],'Interpreter','latex','Fontsize',22)
StorageF=[];
for i = 1:length(sol.value(x))
    StorageF = [StorageF  Storage_Function(sol.value(x(:,i)),sol.value(u(i)))];
end
ffff = figure(5)
plot(time,StorageF)
legend(['$\Psi(x)$'],'Interpreter','latex','Fontsize',26)
xlabel(['Time'])
%% SAVE
XSAVE = sol.value(x);
USAVE = sol.value(u);
DSAVE = sol.value(d);
SFSAVE = StorageF;
save ('100secondsOCPLowParetoPoint.mat', 'XSAVE','USAVE','DSAVE','SFSAVE','time')
%%  This part runs the normal Pareto optimization. The solutions
%% Are stored in the solEE object. 

solEE = nbi2d( ocp, ep_ocp, p_params, nPoints );
 figure(12)
 for i = 1:nPoints
scatter(solEE(i).value(x(6,end)),solEE(i).value(x(7,end)),200,'filled')
hold on
xlabel('Energy')
ylabel('Damage')
end