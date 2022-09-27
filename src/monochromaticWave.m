function waveDisturbance = monochromaticWave(t, T_wave, H_wave)
persistent Gamma_lu

if nargin < 3
    H_wave = 2;
end
if nargin < 2
    T_wave = 10;
end

if isempty(Gamma_lu)
    load('PolySurge_inputs.mat', 'Gamma_lu');
end

waveDisturbance = @(t) 0.5 * H_wave *interp1(Gamma_lu.omega, Gamma_lu.Gamma, 2 * pi / T_wave) * sin(2 * pi * t / T_wave);
end