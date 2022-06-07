%% MCP WITH ClOSE LOOK AT THE DIFFERENT SOLUTION
%% START WITH NORMAL NBI TO FIND GOOD PARAMETERS FOR COST FUNCTION
timehorizon = 50;
timestep = 0.4; 
nPoints = 45;           %% Should be a fast algorythm so you can have lots of points
nCores = 6;             %% Please pay close attention to RAM. If full reduce 

solA = ParalellNBI(timehorizon,timestep,nPoints,nCores) ;

delete(gcp('nocreate')) %% shutdown pool 

%% Plot results for convinient choice of PP
UU = zeros(length(solA),length(solA(1).u));
XX = zeros(length(solA),length(solA(1).x(2,:)));
PP= zeros(length(solA),2);
figure(1)

for i = 1:length(solA)
UU(i,:)=solA(i).u;
XX(i,:)=solA(i).x(2,:);
PP(i,:)=[-solA(i).x(6,end) solA(i).x(7,end)];
% disp([num2str(solA(i).time),solA(i).warmstart])
scatter(PP(i,1),PP(i,2),100,'filled');
hold on
 text(PP(i,1), PP(i,2),['\leftarrow',int2str(i)],'Interpreter','tex');
end
title('Pareto Front')
xlabel(['Energy'])
ylabel('damage')

%% Please chose a point from the pareto Front and then inspect it. 
Pindx = 16; %% Manually change this!
time = [0:timestep:timehorizon];
figure(12)
subplot(1,2,1)
plot(time,UU(Pindx,:))
xlabel('time')
ylabel('Voltage')

subplot(1,2,2)
plot(time,XX(Pindx,:))
xlabel('time')
ylabel('angle')
%% Start MCP here
timehorizon = 50;
timestep = 0.4;
[ocp,x,u,d,x0] = initializeOCP_NBI(timehorizon,timestep,solA(Pindx).ParetoParameters);

sol=ocp.solve();
uOld1= sol.value(u);
xOld1= sol.value(x);
dOld1= sol.value(d);
%% Start MPC
tic;
timeshift = 7;
timeshiftindx= round(length(x)/(timehorizon/timeshift));
x0new =xOld(:,timeshiftindx);
timeshift1 = time(timeshiftindx);
initialGuessX = [xOld1(:,timeshiftindx:end) repmat(xOld1(:,end),1,timeshiftindx-1)];
initialGuessU = [uOld1(timeshiftindx:end) repmat(uOld1(end),1,timeshiftindx-1)];


xOld = xOld1;
uOld = uOld1;

timeold=time;
for i = 1:20
timeshift = timeshift1*i;
[ocp,x,u,d,x0] = initializeOCP_NBI(timehorizon,timestep,solA(Pindx).ParetoParameters,x0=x0new,ts=timeshift);
ocp.set_value(d,arrayfun(@(t) StochasticWave(t),[timestep+timeshift:timestep:(d.length()*timestep)+timeshift]));
ocp.set_initial(x,initialGuessX);
ocp.set_initial(u,initialGuessU);
sol1 = ocp.solve();
time1 =[timestep+timeshift:timestep:(d.length()*timestep)+timeshift];
f = figure(1);
subplot(2,1,1)
plot(timeold,uOld,"Color",'r','LineStyle','-','LineWidth',5)
hold on 
 plot(time1,sol1.value(u),"Color",'b','LineStyle','--','LineWidth',5)
hold off
ylim([0 1200]);
xlim([timeold(1) time1(end)]);
subplot(2,1,2)
plot(timeold,xOld(2,:),"Color",'r','LineStyle','-','LineWidth',5)
hold on 
plot(time1,sol1.value(x(2,:)),"Color",'b','LineStyle','--','LineWidth',5)
hold off
ylim([-1.5 1.5]);
xlim([timeold(1) time1(end)]);
saveU = getframe(f).cdata;
cd Images/
filename =[num2str(i),'_U.png'];
imwrite(saveU ,filename);
cd ..
% ff = figure(2);
% 
% plot(timeold,xOld(2,:),"Color",'r','LineStyle','-','LineWidth',5)
% hold on 
% plot(time1,sol1.value(x(2,:)),"Color",'b','LineStyle','--','LineWidth',5)
% hold off
% saveX = getframe(ff).cdata;
% filename =[num2str(i),'_X.png'];
% imwrite(saveX ,filename);
xOld = sol1.value(x);
uOld = sol1.value(u);
x0new =xOld(:,timeshiftindx);
timeold=time1;
initialGuessX = [sol1.value(x(:,timeshiftindx:end)) repmat(sol1.value(x(:,end)),1,timeshiftindx-1)];
initialGuessU = [sol1.value(u(timeshiftindx:end)) repmat(sol1.value(u(end)),1,timeshiftindx-1)];

end
disp('total time with warmstart was')
disp(num2str(toc))