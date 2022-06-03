function waveDisturbance = StochasticWave(t, T_wave, H_wave,N_freq,d_w)
persistent Gamma_lu
global seedPhrase
global Gamma_nu
global IndividualWaves
seedPhrase = 1;

rng(seedPhrase,'philox');

if nargin < 5
    d_W = pi/20;  %Discretization of Omega
end

if nargin < 4
    N_freq = 10;  %Number of Frequencies in Wave
end

if nargin < 3
    H_wave = 90;  % Hight factor

end
if nargin < 2
    T_wave = 10; % Period of dominant Wave 
end
w_end = 2*2*pi*(1/T_wave);
% w_end = 0.99999;
if isempty(Gamma_lu)
    load('PolySurge_inputs.mat', 'Gamma_nu');
end
waveSpectrum = @(w) 262.9*H_wave^2*T_wave^(-4)*w^(-5)*exp(-1054*T_wave^(-4)*w^(-4));

DataSample = (0.05:d_W:w_end);
Ak =@(w) (2*d_W*waveSpectrum(w))^0.5;
% Gamma = @ (w)  interp1(Gamma_lu.omega, Gamma_lu.Gamma, w);
Gamma_nu = [[0  Gamma_nu(1,2)] ;Gamma_nu];
Gamma = @(w) interp1(Gamma_nu(:,1), Gamma_nu(:,2), w);
% randomShift = datasample(0:pi/40:2*pi,length(DataSample));
randomShift = rand(1,length(DataSample))*2*pi;
% waveDisturbance = sum(arrayfun(Ak,DataSample).*arrayfun(Gamma,DataSample).*sin(DataSample.*t+randomShift));
waveDisturbance=0;
IndividualWaves = zeros(length(DataSample),1);
for i=1:length(DataSample)
    waveDisturbance = waveDisturbance+Ak(DataSample(i)).*Gamma(DataSample(i)).*sin(DataSample(i).*t+randomShift(i));
    IndividualWaves(i)= Ak(DataSample(i)).*Gamma(DataSample(i)).*sin(DataSample(i).*t+randomShift(i));
end


end