% Parameters
Hs = 3;           % Significant wave height (meters)
Tp = 5;           % Peak period (seconds)
gamma = 3.3;      % Peakedness parameter

% Frequency range (based on the peak period)
f_min = 0.05;     % Minimum frequency (Hz)
f_max = 1/Tp;     % Maximum frequency (based on the peak period)

% Simulation parameters
T_total = 600;    % Total simulation time (seconds)
dt = 0.1;         % Time step (seconds)
t = 0:dt:T_total; % Time vector

% JONSWAP Spectrum Parameters
omega_p = 2*pi/Tp;     % Peak angular frequency (rad/s)
sigma = zeros(size(t));

% Frequency array (angular frequencies)
frequencies = linspace(f_min, f_max, length(t));
angular_freq = 2 * pi * frequencies;

% JONSWAP spectrum function
S = @(f) (Hs^2 / 16) * (gamma^exp(-((f - omega_p).^2) / (2 * sigma.^2))) ...
         .* exp(-1.25 * (f / omega_p)^(-4));

% Generate random phases
phi = 2 * pi * rand(1, length(angular_freq));

% Generate wave elevation (η(t)) based on the JONSWAP spectrum
eta = zeros(size(t));
for i = 1:length(frequencies)
    eta = eta + sqrt(2 * S(frequencies(i)) * (frequencies(i))) ...
            * cos(angular_freq(i) * t + phi(i));
end

% Plot wave elevation profile
figure;
plot(t, eta);
title('Wave Elevation Profile');
xlabel('Time (seconds)');
ylabel('Wave Elevation (meters)');
grid on;
