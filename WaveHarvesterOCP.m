global PathToParameters
PathToParameters= 'src/PolySurge_inputs.mat';
load(PathToParameters);




tic;
timehorizon     = 100;         % [1-inf]  How long
timestep        = 0.5;         % [0.05-1] MPC timestep, i.e. the discretisation of the ocp.
                               %          A step larger then 1 is not recommended.
SwingInTime     = 200;         % [100:~]  How long the system is left alone to swing in 
                               %          before the optimal control is applied. 
MOOPMethod      = 'nbi';       % ['nbi' 'ws'] either weighted sum or normal boundary intersection
WaveForm        = 'Stochastic';  % ['Harmonic' 'Stochastic'] Choose the Wave Distrubance
Seed            = 2  ;         % [1-10]   Seed of the Wave distrurbance. 

saving          = true;        % If saving, the results will be saved to the "Results" folder.
plotting        = true;        % If plotting, the resulting paretofront will be plotted

filenameMOOP = ['MOOP' WaveForm '_' num2str(timehorizon) 'seconds_' MOOPMethod '.mat'];             
                               % If saving use this filename

nSteps          = round(timehorizon/timestep);     % Number of discrete timesteps

%create OCP object and apply wave harvester DGL
[ocp,x,u,d,x0_p] = initializeOCPENERGY(timehorizon,timestep); 


time            = linspace(0,timehorizon,d.length());% Create array with discrete time steps
WaveTime        = time+SwingInTime;                  % To create a smooth transition from the swing in the wave 
                                                     % continous after the swing in
                                                      

% Swing in the system for x seconds and set the initial value
x0 = SwingIn(SwingInTime,WaveForm,Seed=Seed);
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
costfun = ([x(6,end) x(7,end)]);

switch MOOPMethod
    case 'ws'
        ptDistance = 0.05;
        [p_params, ep, norm_costfun] = scalarize_moocp( ocp, costfun, method="ws" );
        [solEE, costs] = awds( ocp, norm_costfun, ep, p_params ,ptDistance);
       

    case 'nbi'
        nPoints = 20;
        [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method="ps", normalize='fix' );
        solEE = nbi2d( ocp, ep_ocp, p_params, nPoints );

end
sol = cell(length(solEE),1);
for i = 1:length(solEE)
    sol{i}.x = solEE(i).value(x);
    sol{i}.u = solEE(i).value(u);
    sol{i}.time = time;
    sol{i}.d = solEE(i).value(d);
end

if (plotting)
    f = figure(1)
  
    for i = 1:length(sol)
        scatter(sol{i}.x(6,end),sol{i}.x(7,end),250,'filled')
        hold on
    end

    l = legend();
    for i = 1:length(sol)
        l.String{i} = num2str(i);

    end
    xlabel('Energy')
    ylabel('Damage')
    title (['Pareto Front using ' MOOPMethod ' Method']);
    EGFixFigure
end
if (saving)
    if ~exist([pwd filesep 'Results'],'dir')
       mkdir('Results')
    end
    save(['Results' filesep filenameMOOP],"sol")
end



