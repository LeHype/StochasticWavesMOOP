function [waveDisturbance varargout]= NewStochasticWave(t, args)


arguments    
    t               (1,1) {mustBeNumeric}
    args.T_wave     (1,1) {mustBeNumeric} = 10  % Period of dominant Wave 
    args.H_wave     (1,1) {mustBeNumeric}= 90   % Hight factor
    args.N_freq     (1,1) {mustBeNumeric}= 10   %Number of Frequencies in Wave
    args.d_W        (1,1) {mustBeNumeric}= pi/20 
    args.randomseed (1,1) {mustBeNumeric} = 1
    args.DrawDebug  (1,1) {boolean} = false;
    
end
if args.DrawDebug
    increment = 0.5;
    f = figure(1)

    x1= [0 1];
    y1= [NewStochasticWave(0) NewStochasticWave(increment)];
    p = plot(x1,y1,'LineWidth',6)

    x2 = [];
    y2 = [];

p.XDataSource = 'x1';
p.YDataSource = 'y1';
k = 2*increment;
while(args.DrawDebug)
    [wave Individual]= NewStochasticWave(k);
    x1 = [x1 k];
    y1 = [y1 wave];
    k = k+increment;
    x2 = [x2 k];
    y2 = [y2 Individual];
    refreshdata(f,'caller')
    drawnow

%     pause(0.2)
end
end

persistent Gamma_lu

seedPhrase = args.randomseed;

rng(seedPhrase,'philox');

d_W = args.d_W;
N_freq = args.N_freq;  
H_wave = args.H_wave;
T_wave = args.T_wave;




w_max= 2*2*pi*(1/T_wave);  %% What Range of frequencies will be present in Wave
w_min= 0.05;

if isempty(Gamma_lu)
    load('PolySurge_inputs.mat', 'Gamma_nu'); 
end
waveSpectrum = @(w) 262.9*H_wave^2*T_wave^(-4)*w^(-5)*exp(-1054*T_wave^(-4)*w^(-4));

DataSample = (w_min:(w_max-w_min)/(N_freq-1):w_max);
Ak =@(w) (2*d_W*waveSpectrum(w))^0.5;

Gamma_nu = [[0  Gamma_nu(1,2)] ;Gamma_nu];
Gamma = @(w) interp1(Gamma_nu(:,1), Gamma_nu(:,2), w);


randomShift = rand(1,length(DataSample))*2*pi;


waveDisturbance=0;
IndividualWaves = zeros(length(DataSample),1);
for i=1:length(DataSample)
    waveDisturbance = waveDisturbance+Ak(DataSample(i)).*Gamma(DataSample(i)).*sin(DataSample(i).*t+randomShift(i));
    IndividualWaves(i)= Ak(DataSample(i)).*Gamma(DataSample(i)).*sin(DataSample(i).*t+randomShift(i));

end
varargout= cell(1,1);
varargout{1} = IndividualWaves;

end