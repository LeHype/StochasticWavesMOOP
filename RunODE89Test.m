function [x_hist,u,d_out,du,t] = RunODE89Test(sol,params)
% addpath("src\")
% clearvars


u_dt = sol.time(2);
x_dt = sol.time(2);

timing = 0:x_dt:sol.time(end);
get_energy = false;


model = get_model(timing, sol.u,params);
Storage_Function = @(x,u) 0.5*Parameters.Mh*x(1)^2 + 0.5*Parameters.Kh*x(2)^2 + 0.5*x(3:5)'*Parameters.Q*x(3:5) + 0.5*(Parameters.C0-Parameters.gamma*x(2)^2)*u; 

[t,x_hist] = ode89(model, [0:sol.time(2):sol.time(end)], sol.x(:,1),odeset('Stats','on')); % ,'OutputFcn',@odeplot
x_hist = x_hist';
% sol.u = sqrt(sol.u);
% sol.u = u = repelem(sol.u,10)
pp = spline(timing,sol.u);
u_der=fnder(pp,1);
du=ppval(u_der,t);
u = ppval(pp,t);
wave = monochromaticWave();
d_out = arrayfun(@(a) wave(a),t)';
    function model = get_model(t_u, u_val,params)
        model = @wecdeg;


        function xdot  = wecdeg(t,x)
            persistent Ac Bc R0 gamma
            if isempty(Ac) || isempty(Bc) || isempty(R0) || isempty(gamma)
                load("PolySurge_inputs" ,'Ac', 'Bc', 'R0', 'gamma');
            end
            d = monochromaticWave();
            pp = spline(t_u,u_val);
            u_der=fnder(pp,1);
            du=ppval(u_der,t);
            u = ppval(pp,t);
            xdot = [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d(t);
                cost_energy(x,u,d(t),du,params);
                d(t)*x(1)*1e-6
                ];
        end




    end
end

