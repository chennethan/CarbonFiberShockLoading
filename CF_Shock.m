clear; clc;

%% ---------------- Material properties (psi) ----------------
mat_properties = [23900000, 23900000, 20000000;  % Longitudinal tensile modulus E1
                   1240000,  1240000,  1200000;  % Transverse tensile modulus E2
                     0.300,    0.326,    0.250;  % Poisson's ratio nu12
                    640000,   640000,   800000;  % Shear modulus G12
                    415000,   463000,   300000;  % Longitudinal tensile strength Xt
                   -209000,  -209000,  -150000;  % Longitudinal compressive strength Xc
                     11900,    11900,     7000;  % Transverse tensile strength Yt
                   -100000,  -100000,   -25000;  % Transverse compressive strength Yc
                     20500,    20500,    14000]; % Shear strength S

material = 1;

E_1   = mat_properties(1, material);
E_2   = mat_properties(2, material);
nu_12 = mat_properties(3, material);
nu_21 = E_2 / E_1 * nu_12;
G_12  = mat_properties(4, material);
X_t   = mat_properties(5, material);
X_c   = mat_properties(6, material);
Y_t   = mat_properties(7, material);
Y_c   = mat_properties(8, material);
S     = mat_properties(9, material);

%% ---------------- Tank geometry & pressure ----------------
t_total = 0.069;        % Thickness [in]
R       = 7.764/2;      % Radius [in]
p       = 750;          % Internal pressure [psi]
l       = 6 * 39.3701;  % Total tank length [in] (not used below)

%% ---------------- Parachute opening (tau_inflate governs) ----------------
S = load('tau_outputs.mat'); % Loading tau_inflate results

% Inputs at deployment (SI here; converted where needed)
m_kg        = 34.36;    % Vehicle mass [kg] (include everything the lines decelerate)
V0          = S.DeltaV;     % Downward speed at opening [m/s] derived by tau_inflate
rho         = 1.0;      % Air density at opening altitude [kg/m^3]
Cd          = 1.5;      % Canopy Cd (full inflation) [-]
A_full      = 2.0;      % Canopy area [m^2]
tau_inflate = S.tau_inflate;      % Half-sine pulse length tau [s] (for reporting/tests)

g = 9.80665;            % [m/s^2]

% Peak drag deceleration (no gravity); then include weight for line force
a_peak_tau = S.a_peak_tau; % peak decel derived by tau_inflate
F_peak_N  = m_kg * (a_peak_tau + g); % [N]
F_peak_lb = F_peak_N * 0.224809; % [lb]

% Cross checking using opening shock factor as sanity check (0.9-1.2 is a
% reasonable target band)
K_osf_implied = 2*(F_peak_N - m_kg*g) / (rho * V0^2 * Cd * A_full);   % [m/s^2]

% Half-sine duration metrics (for documentation/dynamics; not used for membrane stress)
tau       = tau_inflate;                 % [s]
dt_rect   = (2/pi) * tau;                % Rectangular-equivalent [s]
dt_05_95  = 0.713 * tau;                 % 5–95% impulse duration [s]
dt_fwhm   = (2/3) * tau;                 % Full width at half max [s]

% Convert peak line force to axial membrane resultant N_x_shock [lb/in]
C_in     = 2*pi*R;                       % Circumference [in]
Nx_shock = F_peak_lb / C_in;             % [lb/in]  (AXIAL ONLY)

% Thin-cylinder pressure resultants [lb/in]
Nx_pressure = p * R / 2;                 % Axial from pressure
Ny_pressure = p * R;                     % Hoop from pressure

% Total in-plane resultants for CLT
Nx  = Nx_pressure + Nx_shock;            % Chute adds to axial
Ny  = Ny_pressure;                       % Hoop unchanged by chute
Nxy = 0;

N = [Nx; Ny; Nxy];
stresses_global = N / t_total;           % [psi]

%% ---------------- Laminate ABD (membrane only here) ----------------
Q = zeros(3);
Q_den  = 1 - nu_12 * nu_21;
Q(1,1) = E_1 / Q_den;
Q(1,2) = nu_21 * E_1 / Q_den;
Q(2,1) = nu_12 * E_2 / Q_den;
Q(2,2) = E_2 / Q_den;
Q(3,3) = G_12;

% Example sublaminate A-matrices (you can swap which stack to analyze)
A_55 = layer_A(55.0, t_total, Q);
A_45 = layer_A(45.0, 0.045,   Q);
A_86 = layer_A(86.0, 0.006,   Q);
A_15 = layer_A(15.0, 0.018,   Q);

A = A_55;                 % Use your chosen stack here
A_inv = inv(A);

% Laminate effective properties (for reference)
E_x  = 1 / (t_total * A_inv(1,1));
E_y  = 1 / (t_total * A_inv(2,2));
nu_xy = -A_inv(1,2) / A_inv(1,1);
G_xy  = 1 / (t_total * A_inv(3,3));

% Membrane strains
strains_global = A \ N;

%% ---------------- Tsai–Wu at +55/-55 plies ----------------
[TW_55,    FS_55,    G_stresses_55,    L_stress_pos55]  = layer_TW( 55.0,  strains_global, Q, mat_properties, material);
[TW_neg55, FS_neg55, G_stresses_neg55, L_stress_neg55]  = layer_TW(-55.0,  strains_global, Q, mat_properties, material);

%% ---------------- Reporting ----------------
disp('--- Parachute Opening (tau_inflate governs) ---');
fprintf('a_peak_tau = %.3f m/s^2\n', a_peak_tau);
fprintf('F_peak = %.0f N (%.0f lb)\n', F_peak_N, F_peak_lb);
fprintf('Durations: tau = %.3f s | dt_rect = %.3f s | dt_05-95 = %.3f s | dt_fwhm = %.3f s\n', ...
        tau, dt_rect, dt_05_95, dt_fwhm);
fprintf('Calculated opening shock factor (Implied) = %.3f', K_osf_implied);

fprintf('\nMembrane resultants [lb/in]: N_x_shock = %.2f | N_x_total = %.2f | N_y_total = %.2f\n', ...
        Nx_shock, Nx, Ny);
fprintf('Global membrane stresses [kpsi]: Axial = %.3f | Hoop = %.3f | Shear = %.3f\n', ...
        stresses_global(1)/1000, stresses_global(2)/1000, stresses_global(3)/1000);

disp('--- Local Ply Stresses [kpsi] ---');
disp('+55° lamina');
fprintf('  Longitudinal: %.3f  |  Transverse: %.3f  |  Shear: %.3f\n', ...
        L_stress_pos55(1)/1000, L_stress_pos55(2)/1000, L_stress_pos55(3)/1000);
disp('-55° lamina');
fprintf('  Longitudinal: %.3f  |  Transverse: %.3f  |  Shear: %.3f\n', ...
        L_stress_neg55(1)/1000, L_stress_neg55(2)/1000, L_stress_neg55(3)/1000);

disp('--- Tsai–Wu ---');
fprintf('+55°: TW_total=%.3f  |  FS=%.3f (reserve factor)\n', TW_55(1), FS_55);
fprintf('-55°: TW_total=%.3f  |  FS=%.3f (reserve factor)\n', TW_neg55(1), FS_neg55);

%% ---------------- Helper functions ----------------
function A = layer_A(wind_angle, thickness, Q)
    theta_pos = deg2rad(wind_angle);
    theta_neg = -theta_pos;

    m_pos = cos(theta_pos); n_pos = sin(theta_pos);
    T_1_pos = [     m_pos^2     n_pos^2   2*m_pos*n_pos;
                    n_pos^2     m_pos^2  -2*m_pos*n_pos;
               -m_pos*n_pos m_pos*n_pos m_pos^2-n_pos^2];
    T_2_pos = [       m_pos^2       n_pos^2     m_pos*n_pos;
                      n_pos^2       m_pos^2    -m_pos*n_pos;
               -2*m_pos*n_pos 2*m_pos*n_pos m_pos^2-n_pos^2];
    Q_bar_pos = T_1_pos \ Q * T_2_pos;

    m_neg = cos(theta_neg); n_neg = sin(theta_neg);
    T_1_neg = [     m_neg^2     n_neg^2   2*m_neg*n_neg;
                    n_neg^2     m_neg^2  -2*m_neg*n_neg;
               -m_neg*n_neg m_neg*n_neg m_neg^2-n_neg^2];
    T_2_neg = [       m_neg^2       n_neg^2     m_neg*n_neg;
                      n_neg^2       m_neg^2    -m_neg*n_neg;
               -2*m_neg*n_neg 2*m_neg*n_neg m_neg^2-n_neg^2];
    Q_bar_neg = T_1_neg \ Q * T_2_neg;

    A = Q_bar_pos * thickness / 2 + Q_bar_neg * thickness / 2;
end

function [TW, FS, G_stresses_pos, stresses_pos] = layer_TW(wind_angle, strains_global, Q, mat_properties, material)
    theta_pos = deg2rad(wind_angle);
    m_pos = cos(theta_pos); n_pos = sin(theta_pos);

    T_1_pos = [     m_pos^2     n_pos^2   2*m_pos*n_pos;
                    n_pos^2     m_pos^2  -2*m_pos*n_pos;
               -m_pos*n_pos m_pos*n_pos m_pos^2-n_pos^2];
    T_2_pos = [       m_pos^2       n_pos^2     m_pos*n_pos;
                      n_pos^2       m_pos^2    -m_pos*n_pos;
               -2*m_pos*n_pos 2*m_pos*n_pos m_pos^2-n_pos^2];

    strains_pos  = T_2_pos * strains_global;
    stresses_pos = Q * strains_pos;        % local (1-2-12) stresses
    sigma_1 = stresses_pos(1,1);
    sigma_2 = stresses_pos(2,1);
    tau_12  = stresses_pos(3,1);

    G_stresses_pos = T_1_pos \ stresses_pos; % optional: back to global

    X_t = mat_properties(5, material);
    X_c = mat_properties(6, material);
    Y_t = mat_properties(7, material);
    Y_c = mat_properties(8, material);
    S   = mat_properties(9, material);

    % Tsai–Wu coefficients
    F_11 = -1 / (X_t * X_c);
    F_1  =  1 / X_t + 1 / X_c;
    F_22 = -1 / (Y_t * Y_c);
    F_2  =  1 / Y_t + 1 / Y_c;
    F_66 =  1 / S^2;

    % Tsai–Wu index components
    TW_total = F_1*sigma_1 + F_2*sigma_2 ...
             + F_11*sigma_1^2 + F_22*sigma_2^2 + F_66*tau_12^2 ...
             - F_11*sigma_1*sigma_2;
    TW_fiber = F_1*sigma_1 + F_11*sigma_1^2;
    TW_resin = F_2*sigma_2 + F_22*sigma_2^2 + F_66*tau_12^2 - F_11*sigma_1*sigma_2;

    TW = [TW_total; TW_fiber; TW_resin];

    % Reserve factor (load multiplier to failure) from TW=1 quadratic
    FS_A = F_11*sigma_1^2 + F_22*sigma_2^2 + F_66*tau_12^2 - F_11*sigma_1*sigma_2;
    FS_B = F_1*sigma_1 + F_2*sigma_2;
    FS   = (sqrt(FS_B^2 + 4*FS_A) - FS_B) / (2*FS_A);
end