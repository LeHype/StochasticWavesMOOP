function [solA] = MPCNBI(timehorizon,timestep,nPoints,Wave,args)
arguments
               
    timehorizon  (1,1) {mustBeNumeric}
    timestep     (1,1) {mustBeNumeric}
    nPoints      (1,1) {mustBeNumeric}
    Wave         (:,:) {mustBeNumeric}
    args.x0      (7,1) {mustBeNumeric} = zeros(7,1)
    args.ts      (1,1) {mustBeNumeric} = 0.0
    args.u1      (1,1) {mustBeNumeric} = 0.0
   end
[ocp,x,u,d,x0,x0_p] = initializeOCP(timehorizon,timestep, x0 = args.x0);
ocp.set_value(u(1),args.u1);
ocp.set_value(x0_p,x0);
ocp.set_value(d,Wave);
NumInc =round(timehorizon/(timestep));

  
costs= ([x(6,end) x(7,end)]);
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


%% Number of iterations per loop 
%This is to ensure warmstarts for more the the initial 8 points



weights = linspace(1,0,nPoints);
weights = [weights; 1-weights];

% bp_pts = weights(1,:)'.*ep(1,:) + weights(2,:)'.*ep(2,:);


solA = [];
startpara=tic;

sdir = (ocp.parameter(2));
l   = ocp.variable( 1 );
pt  = ocp.parameter( 1, 2 );
dir = ocp.parameter( 1, 2 );
for j = 1:nPoints
startinner=tic;


bp_pts = weights(1,:)'.*ep(1,:) + weights(2,:)'.*ep(2,:);

costs= ([x(6,end) x(7,end)]-[up])./[np-up];


ocp.set_value(sdir, ndir);


ocp.subject_to( pt + l*dir >= costs )
        
ocp.minimize(-l)
ocp.set_value(pt,bp_pts(j,:));
ocp.set_value(dir,ndir);
warmstart='false';
if j>1
%       ocp.set_initial([x(:); u(:); ocp.lam_g], solprev.value([x(:); u(:); ocp.lam_g]))
      ocp.set_initial(x(:,2:end),sol.x(:,1:end-1))
      ocp.set_initial(u(2:end),sol.u(2:end))     
%       ocp.set_initial(ocp.lam_g,lamG)
     solprev.delete()
     warmstart='true';
end
 solprev = ocp.solve();
 lamG = ocp.value(ocp.lam_g);
 sol = struct;
 ParetoParameters = struct;
 ParetoParameters.searchdirection=ndir;
 ParetoParameters.startingpoint=bp_pts(j,:);
 ParetoParameters.up = up;
 ParetoParameters.np = np;
 sol.startingpoint=bp_pts(j,:);
 sol.x = solprev.value(x);
 sol.u = solprev.value(u);
 sol.stats = solprev.stats;
 sol.warmstart=warmstart;
 sol.time=toc(startinner);
 sol.ParetoParameters=ParetoParameters;
 solA = [solA sol];
 
 

end


end



