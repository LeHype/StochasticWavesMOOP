function [x0,varargout] = SwingIn(SwingInTime,WaveForm, x0_p, args)
arguments
    SwingInTime  (1,1) {mustBeNumeric}
    WaveForm     (:,:)   {mustBeText}
    x0_p
    args.dt      (1,1) {mustBeNumeric} = 0.5
    args.Seed    (1,1) {mustBeNumeric} = 1
end
global PathToParameters;
load(PathToParameters ,'Ac', 'Bc', 'gamma','R0');
wave_dgl = @(x,u,d) [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d];
                              
x0 = zeros(5,1);
StateHist = zeros(5, SwingInTime/args.dt);
WaveHist = zeros(SwingInTime/args.dt,1);
for i = 1:round(SwingInTime/args.dt)
    switch WaveForm
        case 'Harmonic'
            HarmonicWave = monochromaticWave();
            x0 = integrator_step_disturbed(x0,[0],args.dt,wave_dgl,HarmonicWave(i*args.dt));
            WaveHist(i) = HarmonicWave(i*args.dt);
        case 'Stochastic'
            x0 = integrator_step_disturbed(x0,[0],args.dt,wave_dgl,FBMStochasticWave(i*args.dt,Seed=args.Seed));
            WaveHist(i) = FBMStochasticWave(i*args.dt,Seed=args.Seed);
    end
    StateHist(:,i)=full(evalf(x0(:)));
    
end
x0= [full(evalf(x0)); zeros(size(x0_p,1)-size(x0,1), 1)];
varargout{1} = StateHist;
varargout{2} = WaveHist;

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
end

