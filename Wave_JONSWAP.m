function elevation=Wave_JONSWAP(Height,T,Gamma)
% Parameters for the JONSWAP spectrum
g = 9.81;        % Gravitational acceleration (m/s^2)
H_s = Height;         % Significant wave height (m)
T_z = T;         % Zero-crossing period (s)
alpha = 0.0081;  % JONSWAP spectrum constant for H_s
gamma = Gamma;     % Peak enhancement factor (default 3.3)
sigma = 0.07;    % Spectral width
omega_p = 2*pi/T_z;  % Peak angular frequency

% Define frequency range for the spectrum
f_min = 0.01;    % Minimum frequency (Hz)
f_max = 2;       % Maximum frequency (Hz)
N = 1000;        % Number of frequency bins
df = (f_max - f_min) / N;  % Frequency resolution
f = f_min:df:f_max;  % Frequency vector

% JONSWAP Spectrum calculation
omega = 2 * pi * f;    % Angular frequency
S_omega = alpha * g^2 ./ omega.^5 .* exp(-1.25 * (omega_p ./ omega).^4) .* gamma.^(exp(-((omega - omega_p).^2) / (2*sigma^2)));

% Generate random phases for each frequency bin
phi = 2*pi*rand(1, length(omega));  % Ensure phi has the same length as omega

% Generate complex amplitudes for each frequency bin
A = sqrt(S_omega) .* exp(1i * phi);

% Create time vector
T = 600;  % Total simulation time (seconds)
dt = 0.1; % Time step (seconds)
t = 0:dt:T;  % Time vector

% Calculate wave elevation as the sum of harmonic components
elevation = zeros(size(t));  % Initialize wave elevation
for k = 1:length(omega)
    elevation = elevation + real(A(k) * exp(1i * omega(k) * t));
end

% Normalize the elevation by the significant wave height
elevation = elevation / max(abs(elevation)) * H_s;
end
