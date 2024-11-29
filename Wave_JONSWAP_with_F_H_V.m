% Wave Height, Speed, and Force Calculation
function [wave_elevation,wave_speed,wave_force]=Wave_JONSWAP_with_F_H_V(Height,T,Gamma)
% Parameters
rho = 1025; % Water density (kg/m^3)
g = 9.81; % Gravity (m/s^2)
H_s = Height;%2.0; % Significant wave height (m)
T_p = T;%8.0; % Peak wave period (s)
gamma =Gamma; %3.3; % JONSWAP peak enhancement factor
sigma_a = 0.07; % Sigma below peak
sigma_b = 0.09; % Sigma above peak
depth = 50; % Water depth (m)
A = 1.0; % Assumed cross-sectional area for force calculation (m^2)

% Frequency domain setup
N = 2048; % Number of frequency points (higher for accuracy)
omega_peak = 2 * pi / T_p; % Peak angular frequency (rad/s)
omega = linspace(0.01, 2 * omega_peak, N); % Angular frequency array (rad/s)
delta_omega = omega(2) - omega(1); % Frequency resolution

% JONSWAP Spectrum Function
jonswap_spectrum = @(omega, omega_peak, H_s, gamma) ...
    (5/16) * (H_s^2 / T_p) .* (omega ./ omega_peak).^(-5) ...
    .* exp(-1.25 * (omega ./ omega_peak).^(-4)) ...
    .* gamma.^(exp(-0.5 * ((omega ./ omega_peak - 1) ./ ...
    (omega <= omega_peak) .* sigma_a + (omega > omega_peak) .* sigma_b).^2));

% Calculate JONSWAP Spectrum
S = jonswap_spectrum(omega, omega_peak, H_s, gamma);

% Random phase for wave synthesis
rng(42); % For reproducibility
phi = 2 * pi * rand(1, N); % Random phases

% Time domain setup
t = linspace(0, 1000, 10000); % Simulation time (1000s, 10 Hz sampling rate)
wave_elevation = zeros(size(t)); % Preallocate wave elevation

% Wave Elevation Synthesis
for i = 1:length(t)
    wave_elevation(i) = sum(sqrt(2 * S * delta_omega) .* cos(omega * t(i) + phi));
end

% Calculate wave speed (particle velocity)
k = omega.^2 / g; % Wave number using shallow-water approximation
wave_speed = zeros(size(t));
for i = 1:length(t)
    wave_speed(i) = sum(sqrt(2 * S * delta_omega) .* omega .* sinh(k * depth) .* cos(omega * t(i) + phi));
end

% Calculate wave force
wave_force = rho * g * A .* wave_elevation; % Hydrostatic approximation

% % Reduce Data Without Using `downsample`
% sampling_interval = 100; % Take every 100th sample
% t_down = t(1:sampling_interval:end);
% wave_elevation_down = wave_elevation(1:sampling_interval:end);
% wave_speed_down = wave_speed(1:sampling_interval:end);
% wave_force_down = wave_force(1:sampling_interval:end);
% 
% % Create Table
% wave_data = table(t_down', wave_elevation_down', wave_speed_down', wave_force_down', ...
%     'VariableNames', {'Time', 'WaveHeight', 'WaveSpeed', 'WaveForce'});

% % Save to CSV
% writetable(wave_data, 'wave_data.csv');
% 
% % Plot Results
% figure;
% 
% % Wave Elevation
% subplot(3, 1, 1);
% plot(t_down, wave_elevation_down, 'b', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Wave Height (m)');
% title('Wave Elevation');
% grid on;
% 
% % Wave Speed
% subplot(3, 1, 2);
% plot(t_down, wave_speed_down, 'r', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Wave Speed (m/s)');
% title('Wave Particle Speed');
% grid on;
% 
% % Wave Force
% subplot(3, 1, 3);
% plot(t_down, wave_force_down, 'g', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Wave Force (N)');
% title('Wave Force');
% grid on;
% 
% % Display first 10 rows of the table
% disp(wave_data(1:10, :));
end