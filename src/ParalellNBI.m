function [solA] = ParalellNBI(timehorizon,timestep,nPoints,NumCores)
NumInc =round(timehorizon/(timestep));

[ocp,x,u,d,x0] = initializeOCP(timehorizon,timestep);

ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));
costs= ([-x(6,end) x(7,end)]);
ocp.minimize( costs*[1 ; 1E-4]);
sol1 = ocp.solve();
ep(1,:) = sol1.value(costs);

ocp.minimize( costs*[1E-4 ; 1]);
sol2= ocp.solve();
ep(2,:) = sol2.value(costs);



up = min( ep, [], 1 );
np = max( ep, [], 1 );


ep = (ep - up)./(np - up);
ndir = null(ep(2:end,:)-ep(1,:));

if all(ndir >= 0)
    ndir = -ndir;
end


if isempty(gcp('nocreate'))
    parpool(NumCores);
    NumWorkers = gcp('nocreate').NumWorkers;
else
    NumWorkers = gcp('nocreate').NumWorkers;
    if NumWorkers ~= NumCores
        warning('The Number of specified cores is not equal to the')
        warning('number of cores in the currently active pool')
    end
end

%% Number of iterations per loop 
%This is to ensure warmstarts for more the the initial 8 points
NInner= ceil(nPoints/NumWorkers);
nPoints = NInner*NumWorkers;


weightsO = linspace(1,0,NumWorkers+1);
weightsO = [weightsO; 1-weightsO];

% bp_pts = weights(1,:)'.*ep(1,:) + weights(2,:)'.*ep(2,:);


solA = [];
startpara=tic;
parfor i = 1:NumWorkers

for j = 1:NInner
startinner=tic;
weights = linspace(weightsO(1,i),weightsO(1,i+1),NInner+1);
weights = [weights; 1-weights];

bp_pts = weights(1,:)'.*ep(1,:) + weights(2,:)'.*ep(2,:);

[ocp,x,u,d] = initializeOCP(timehorizon,timestep,x0=x0);


ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep:timestep:d.length()*timestep]));

costs= ([-x(6,end) x(7,end)]-[up])./[np-up];


sdir = (ocp.parameter(2));
ocp.set_value(sdir, ndir);

l   = ocp.variable( 1 );
pt  = ocp.parameter( 1, 2 );
dir = ocp.parameter( 1, 2 );
ocp.subject_to( pt + l*dir >= costs )
        
ocp.minimize(-l)
ocp.set_value(pt,bp_pts(j,:));
ocp.set_value(dir,ndir);
warmstart='false';
if j>1
     ocp.set_initial(solprev.value_variables)
     solprev.delete()
     warmstart='true';
end
 solprev = ocp.solve();
 sol = struct;
 sol.x = solprev.value(x);
 sol.u = solprev.value(u);
 sol.stats = solprev.stats;
 sol.warmstart=warmstart;
 sol.time=toc(startinner);
 solA = [solA sol];
 
 
%  try 
%  sol = UnJankify(solprev);
% sol.warmstart=warmstart;
% sol.time=toc(startinner)
% 
%  solA= [solA sol];
% 
%  catch
%      sol = struct
%      sol.error = 'true'
%      solA = [solA sol;]
%  end
end


end
totaltime=toc(startpara);
filename= [num2str(timehorizon) ,'_',num2str(timestep),'_', num2str(nPoints),'_',num2str(NumWorkers), '_Para.mat'];
save(filename,'solA','totaltime');
disp('total time parallel ');
disp(totaltime);
end

