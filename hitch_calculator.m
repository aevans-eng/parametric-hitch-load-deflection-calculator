%% HITCH CARRIER LOAD & DEFLECTION ANALYSIS
%  Author: Aaron Evans
%  Date:   January 2026
%  Units:  Imperial (lbs, in, psi)
%
%  Inputs:  Payload weight, material properties, tube geometry, distances
%  Outputs: Factor of safety for welds, torsion, and bending (static & dynamic)
%           Deflection at receiver tip and bike rack COG
%
%  Methods:
%    - Weld shear via weld-as-a-line
%    - Crossbar torsion via Bredt-Batho (thin-wall)
%    - Deflection via superposition (cantilever bending + twist)
%    - Dynamic analysis at 4G pothole load

clc; clear; close all; format compact;

%% 1. INPUT PARAMETERS
% =========================================================================

% -- Loading --
P_load = 180;               % Weight at Center of Gravity (lbs)

% -- Material Properties (A36 Steel) --
E = 29e6;                   % Young's Modulus (psi)
G = 11.2e6;                 % Shear Modulus (psi)
yield_strength = 46000;     % Yield Strength (psi)
allowable_shear = yield_strength / sqrt(3);  % Von Mises shear limit

% -- Weld Properties --
weld_leg = 3/16;            % Weld leg size (in)
electrode_str = 70000;      % Electrode tensile strength (E70 series) (psi)
weld_allowable = electrode_str * 0.3;  % Allowable shear stress

% -- Geometry: Crossbar (2x2 Square Tube) --
b_cross = 2.0;              % Width (in)
h_cross = 2.0;              % Height (in)
t_cross = 0.25;             % Wall thickness (in)
L_cross = 38.0;             % Total length (in)

% -- Geometry: Receiver (2.63x2.63 Square Tube) --
b_rec = 2.63;               % Width (in)
h_rec = 2.63;               % Height (in)
t_rec = 0.25;               % Wall thickness (in)

% -- Geometry: Key Distances --
L_pin_mid = 1.87;           % Hitch pin to midline of crossbar
L_edge_pin = 3.13;          % Receiver end to hitch pin
L_stickout = L_edge_pin + L_pin_mid - 0.5*b_cross;
L_edge_mid = L_pin_mid + L_edge_pin;

% -- Geometry: Rack / Payload Position --
L_rack_x = 29.0;            % Horizontal distance, hitch pin to COG
L_rack_y = 36.0;            % Vertical height, hitch pin to COG
L_rack_x_edge = L_rack_x - L_edge_pin;

%% 2. SECTION PROPERTIES
% =========================================================================

% Weld Geometry
weld_throat = weld_leg * 0.707;
L_weld_dist = L_rack_x_edge - L_edge_pin + L_stickout;

% Receiver Moment of Inertia
b_rec_in = b_rec - 2*t_rec;
h_rec_in = h_rec - 2*t_rec;
Ix_rec = (b_rec*h_rec^3)/12 - (b_rec_in*h_rec_in^3)/12;

% Crossbar Torsion Constants (Bredt-Batho Thin-Wall Theory)
b_m = b_cross - t_cross;    % Median-line width
h_m = h_cross - t_cross;    % Median-line height
Am_cross = b_m * h_m;       % Area enclosed by median line
S_cross = 2 * (b_m + h_m);  % Median-line perimeter
J_cross = 4 * Am_cross^2 * t_cross / S_cross;

%% 3. STATIC STRESS ANALYSIS
% =========================================================================

% --- Weld Shear (weld-as-a-line method) ---
Fmax_weld = (2*b_rec^2) * (weld_throat*weld_allowable) / (b_rec + 2*L_weld_dist);

% --- Crossbar Torsion ---
Tq_cross = P_load * (L_rack_x + L_pin_mid);
shear_stress_cross = (P_load * (L_rack_x + L_pin_mid)) / (2 * t_cross * Am_cross);

% --- Output: Static Results ---
fprintf('\n%s\n', repmat('=',1,50));
fprintf(' STATIC ANALYSIS SUMMARY\n');
fprintf('%s\n', repmat('-',1,50));
fprintf(' %-30s : %8.2f lbs\n', 'Static Load Applied', P_load);
fprintf('\n %-20s | %-12s | %-10s \n', 'COMPONENT', 'F.O.S.', 'STATUS');
fprintf('%s\n', repmat('-',1,50));

% Weld check
fos_weld = Fmax_weld / P_load;
if P_load > Fmax_weld
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-10s\n', 'Weld Shear', fos_weld, status);

% Crossbar check
fos_bar = allowable_shear / shear_stress_cross;
if shear_stress_cross > allowable_shear
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-10s\n', 'Crossbar Torsion', fos_bar, status);

%% 4. DEFLECTION ANALYSIS (Superposition)
% =========================================================================

% --- Crossbar Twist ---
L_effective = (L_cross - b_rec) / 2;
phi_cross = (Tq_cross * 0.5) * (L_effective) / (G * J_cross);

% --- Receiver Bending (Point Load + Moment) ---
delta_rec = (P_load*L_edge_mid^3)/(3*E*Ix_rec) + ...
            (P_load*L_rack_x_edge*L_edge_mid^2)/(2*E*Ix_rec);

% --- Combined Deflection at Receiver Tip ---
delta_total_edge = delta_rec + L_edge_mid * sin(phi_cross);

% --- Tip Angle (Bending Slope + Twist) ---
theta_force_edge = (P_load * L_edge_mid^2) / (2 * E * Ix_rec);
M_at_edge = P_load * L_rack_x_edge;
theta_moment_edge = (M_at_edge * L_edge_mid) / (E * Ix_rec);
phi_total_edge = phi_cross + theta_force_edge + theta_moment_edge;

% --- Deflection at Bike Rack COG ---
delta_bike_rack = delta_total_edge + (L_rack_x_edge * tan(phi_total_edge));

% --- Output: Deflection ---
fprintf('%s\n', repmat('-',1,50));
fprintf(' STATIC DEFLECTION RESULTS\n');
fprintf(' %-30s : %8.2f deg\n', 'Tip Angle (Twist+Bend)', rad2deg(phi_total_edge));
fprintf(' %-30s : %8.4f in\n', 'Vertical Drop (Receiver)', delta_total_edge);
fprintf('      > Due to Bend    : %8.4f in\n', delta_rec);
fprintf('      > Due to Twist   : %8.4f in\n', L_edge_mid * sin(phi_cross));
fprintf(' %-30s : %8.4f in\n', 'Vertical Drop (Bike COG)', delta_bike_rack);
fprintf('%s\n', repmat('=',1,50));

%% 5. DYNAMIC ANALYSIS (Pothole Impact)
% =========================================================================

G_factor = 4.0;
P_dynamic = P_load * G_factor;

fprintf('\n%s\n', repmat('=',1,50));
fprintf(' DYNAMIC ANALYSIS (%.1f G Pothole Load)\n', G_factor);
fprintf('%s\n', repmat('-',1,50));
fprintf(' %-30s : %8.2f lbs\n', 'Dynamic Load', P_dynamic);
fprintf('\n %-20s | %-12s | %-8s \n', 'COMPONENT', 'F.O.S.', 'STATUS');
fprintf('%s\n', repmat('-',1,50));

% Weld (dynamic)
fos_weld_dyn = Fmax_weld / P_dynamic;
if P_dynamic > Fmax_weld
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-8s\n', 'Weld Shear', fos_weld_dyn, status);

% Crossbar torsion (dynamic)
shear_stress_dyn = shear_stress_cross * G_factor;
fos_bar_dyn = allowable_shear / shear_stress_dyn;
if shear_stress_dyn > allowable_shear
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-8s\n', 'Crossbar Torsion', fos_bar_dyn, status);

% Receiver bending stress (dynamic)
M_dynamic_root = P_dynamic * (L_stickout + L_rack_x);
c_dist = h_rec / 2;
sigma_dynamic = (M_dynamic_root * c_dist) / Ix_rec;
fos_rec_dyn = yield_strength / sigma_dynamic;

if sigma_dynamic > yield_strength
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-8s\n', 'Receiver Bending', fos_rec_dyn, status);

fprintf('%s\n', repmat('-',1,50));

% Deflection (dynamic)
delta_dynamic_total_edge = delta_total_edge * G_factor;
delta_dynamic_drop = delta_bike_rack * G_factor;

fprintf(' %-30s : %8.4f in\n', 'Dynamic Tip Drop', delta_dynamic_total_edge);
fprintf('      > Due to Bend    : %8.4f in\n', delta_rec * G_factor);
fprintf('      > Due to Twist   : %8.4f in\n', (L_edge_mid * sin(phi_cross)) * G_factor);
fprintf(' %-30s : %8.4f in\n', 'Dynamic COG Drop', delta_dynamic_drop);
fprintf('%s\n', repmat('=',1,50));
