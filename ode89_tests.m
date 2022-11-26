addpath("src\")
clearvars
load('C:\Work\Projects\StochasticWavesMOOP\src\PolySurge_inputs.mat')

u_dt = 2;
x_dt = 0.1;

timing = 0:x_dt:300;
get_energy = false;
% u_hist = rand(1, length(timing)).^4 * 33000^2;
% u_hist_fine = [repelem(u_hist(1:end-1), u_dt/x_dt), u_hist(end)];
% u_hist_fine = interp1(timing, u_hist, 0:x_dt:300, 'cubic', 'extrap');
u_hist_fine = sin(4*pi*timing);

model = get_model(timing, u_hist_fine, get_energy);
Storage_Function = @(x,u) 0.5*Mh*x(1)^2 + 0.5*Kh*x(2)^2 + 0.5*x(3:5)'*Q*x(3:5) + 0.5*(C0-gamma*x(2)^2)*u; 

[t,x_hist] = ode89(model, 0:x_dt:300, [zeros(5,1); 1e-6*Storage_Function(zeros(5,1),u_hist_fine(1));zeros(get_energy*6,1)]);
x_hist = x_hist';

SF = [];
for i = 1:length(x_hist)
    SF(i) = Storage_Function(x_hist(:,i), u_hist_fine(:,i));
end
plot(diff(SF))
hold on
plot(diff(x_hist(6,:))*1e6)
hold off


function model = get_model(t_u, u_val, with_energy)

if ~with_energy
    model = @wecdeg;
else
    model = @wecdeg_energy;
end

    function xdot = wecdeg(t,x)
        persistent Ac Bc R0 gamma
        if isempty(Ac) || isempty(Bc) || isempty(R0) || isempty(gamma)
            load("PolySurge_inputs" ,'Ac', 'Bc', 'R0', 'gamma');
        end
        d = monochromaticWave();
        u = interp1(t_u, u_val, t, 'linear', 'extrap');
        du = (interp1([t_u, 1e10], [u_val u_val(end)], t+1e-4, 'next', 'extrap') - interp1(t_u, u_val, t, 'previous', 'extrap'))/2;
        xdot = [Ac * x(1:5) - Bc * u * gamma * x(2) + Bc * d(t);
            cost_energy(x,u,du,d(t))*1e-6;
            ];
    end

    function xdot = wecdeg_energy(t,x)
        persistent Ac Bc R0 gamma
        if isempty(Ac) || isempty(Bc) || isempty(R0) || isempty(gamma)
            load("PolySurge_inputs" ,'Ac', 'Bc', 'R0', 'gamma');
        end
        d = monochromaticWave();
        u = interp1(t_u, u_val, t, 'linear', 'extrap');
        du = (interp1([t_u, 1e10], [u_val u_val(end)], t+1e-4, 'next', 'extrap') - interp1(t_u, u_val, t, 'previous', 'extrap'))/0.1;
        xdot = [Ac * x(1:5) - Bc * u * gamma * x(2) + Bc * d(t);
            cost_energy(x, u, du,d(t));
            cost_damage(x,u);
            EE1(x,u,d(t));
            EE2(x,u,d(t));
            EE3(x,u,d(t));
            EE4(x,u,d(t));
            gamma*x(2)*u-u/(R0)
            ];

end
end

function cost = cost_energy(x, u, du, d)
    persistent Ch S R0 C0 gamma
    if isempty(Ch) || isempty(S) || isempty(R0) || isempty(C0) || isempty(gamma)
        load("PolySurge_inputs" ,'Ch', 'S', 'R0', 'C0', 'gamma');
    end
    cost = -(Ch*x(1).^2 + x(3:5)'*S*x(3:5)) + x(1)*(d - 2*gamma.*x(2)*u) + 0.5*(C0-gamma*x(2)^2)*du;
end
function E1 = EE1(x, u, d)
    persistent Ch S R0
    if isempty(Ch) || isempty(S) || isempty(R0)
        load("PolySurge_inputs", 'Ch', 'S', 'R0');
    end
    E1 = 1e-6*(Ch*x(1).^2);
end
function E2 = EE2(x, u, d)
    persistent Ch S R0
    if isempty(Ch) || isempty(S) || isempty(R0)
        load("PolySurge_inputs", 'Ch', 'S', 'R0');
    end
     E2 = (x(3:5)'*S*x(3:5))*1e-6 ;
end
function E3 = EE3(x, u, d)
    persistent Ch S R0
    if isempty(Ch) || isempty(S) || isempty(R0)
        load("PolySurge_inputs", 'Ch', 'S', 'R0');
    end
      E3 = - (d .* x(1))*1e-6;
end
function E4 = EE4(x, u, d)
    persistent Ch S R0
    if isempty(Ch) || isempty(S) || isempty(R0)
        load("PolySurge_inputs", 'Ch', 'S', 'R0');
    end
     E4 = u/R0;
end

function cost = cost_damage(x, u, ~, ~)   
%     cost = (max(u - 484/(cos(x(2)).^2), 0).^2)*1e-6;
    cost = (max(u - 484, 0).^2)*1e-6;
end

