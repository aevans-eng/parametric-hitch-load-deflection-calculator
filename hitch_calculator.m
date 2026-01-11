%% HITCH CARRIER DEFLECTION & STRESS ANALYSIS
%  Author: Aaron Evans
%  Date:   January 2026
%  Units:  Imperial (lbs, in, psi)
%
%  Description: 
%  Analytical verification of a custom fabricated vehicle hitch. 
%  Calculates static stresses on weld, and deflection via 
%  superposition, and dynamic loads based on G-force inputs.

clc; clear; close all; format compact;

%% 1. INPUT VARIABLES

% -- Loading --
P_load = 180;             % Weight at Center of Gravity (lbs)

% -- Material Properties (A36 Steel) --
E = 29e6;                 % Young's Modulus (psi)
G = 11.2e6;               % Shear Modulus (psi)
yield_strength = 46000;   % Yield Strength (psi)
allowable_shear = yield_strength / sqrt(3); % Von Mises shear limit

% -- Weld Properties --
weld_leg = 3/16;          % Weld leg size (in)
electrode_str = 70000;    % Electrode tensile strength (E70 series) (psi)
weld_allowable = electrode_str * 0.3; % Allowable shear stress

% -- Geometry: Crossbar --
b_cross = 2.0;            % Width (in)
h_cross = 2.0;            % Height (in)
t_cross = 0.25;           % Wall thickness (in)
L_cross = 38.0;           % Total length (in)

% -- Geometry: Receiver --
b_rec = 2.63;             % Width (in)
h_rec = 2.63;             % Height (in)
t_rec = 0.25;             % Wall thickness (in)

% -- Geometry: Distances --
L_pin_mid = 1.87;         % Hitch pin to midline of crossbar
L_edge_pin = 3.13;        % Receiver end to hitch pin
L_stickout = L_edge_pin + L_pin_mid - 0.5*b_cross; % Receiver length from end to weld
L_edge_mid = L_pin_mid + L_edge_pin;

% -- Geometry: Rack / Payload --
L_rack_x = 29.0;          % Horizontal distance (Hitch pin to COG)
L_rack_y = 36.0;          % Vertical height (Hitch pin to COG)
L_rack_x_edge = L_rack_x - L_edge_pin;

%% 2. GEOMETRIC PROPERTY CALCULATIONS

% Weld Geometry
weld_throat = weld_leg * 0.707;
L_weld_dist = L_rack_x_edge - L_edge_pin + L_stickout; % Moment arm to weld

% Receiver Inertia
b_rec_in = b_rec - 2*t_rec;
h_rec_in = h_rec - 2*t_rec;
Ix_rec = (b_rec*h_rec^3)/12 - (b_rec_in*h_rec_in^3)/12;

% Crossbar Torsion Constants (Bredt-Batho / Thin-wall)
b_m = b_cross - t_cross;      % Median-line width
h_m = h_cross - t_cross;      % Median-line height
Am_cross = b_m * h_m;         % Area enclosed by median line
S_cross  = 2 * (b_m + h_m);   % Median-line perimeter
J_cross  = 4 * Am_cross^2 * t_cross / S_cross; 

%% 3. STATIC STRESS ANALYSIS

% A. Weld Shear
% Max force on weld group based on treating weld as a line
Fmax_weld = (2*b_rec^2) * (weld_throat*weld_allowable) / (b_rec + 2*L_weld_dist); 

% B. Crossbar Torsion
% Torque applied to crossbar by the rack payload
Tq_cross = P_load * (L_rack_x + L_pin_mid);
% Shear stress in the crossbar walls
shear_stress_cross = (P_load * (L_rack_x + L_pin_mid)) / (2 * t_cross * Am_cross); 

% Output: Static Results
fprintf('\n%s\n', repmat('=',1,50));
fprintf(' STATIC ANALYSIS SUMMARY\n');
fprintf('%s\n', repmat('-',1,50));
fprintf(' %-30s : %8.2f lbs\n', 'Static Load Applied', P_load);
fprintf('\n %-20s | %-12s | %-10s \n', 'COMPONENT', 'F.O.S.', 'STATUS');
fprintf('%s\n', repmat('-',1,50));

% Check 1: Weld
fos_weld = Fmax_weld / P_load;
if P_load > Fmax_weld
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-10s\n', 'Weld Shear', fos_weld, status);

% Check 2: Crossbar
fos_bar = allowable_shear / shear_stress_cross;
if shear_stress_cross > allowable_shear
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-10s\n', 'Crossbar Torsion', fos_bar, status);


%% 4. DEFLECTION CALCULATIONS (Superposition)

% A. Crossbar Twist
% Effective length for torsion (assuming stiffness increase at weld junction)
L_effective = (L_cross - b_rec) / 2; 
% Twist angle (Radians). *0.5 assumes torque is split between left/right arms
phi_cross = (Tq_cross * 0.5) * (L_effective) / (G * J_cross);  

% B. Receiver Bending
% Deflection due to Point Load + Deflection due to Moment
delta_rec = (P_load*L_edge_mid^3)/(3*E*Ix_rec) + ...
            (P_load*L_rack_x_edge*L_edge_mid^2)/(2*E*Ix_rec);

% C. Total Deflection at Receiver Tip
% 1. Vertical Drop (Twist Drop + Bending Drop)
delta_total_edge = delta_rec + L_edge_mid * sin(phi_cross);

% 2. Angle at Tip (Slope due to Force + Slope due to Moment + Twist)
theta_force_edge = (P_load * L_edge_mid^2) / (2 * E * Ix_rec);
M_at_edge = P_load * L_rack_x_edge; 
theta_moment_edge = (M_at_edge * L_edge_mid) / (E * Ix_rec);
phi_total_edge = phi_cross + theta_force_edge + theta_moment_edge;

% D. Total Deflection at Bike Rack COG
% Project linear drop from tip angle
delta_bike_rack = delta_total_edge + (L_rack_x_edge * tan(phi_total_edge));

% Output: Deflection
fprintf('%s\n', repmat('-',1,50));
fprintf(' [3] STATIC DEFLECTION RESULTS\n');
fprintf(' %-30s : %8.2f deg\n', 'Tip Angle (Twist+Bend)', rad2deg(phi_total_edge));
fprintf(' %-30s : %8.4f in\n', 'Vertical Drop (Receiver)', delta_total_edge);
fprintf('      > Due to Bend    : %8.4f in\n', delta_rec);
fprintf('      > Due to Twist   : %8.4f in\n', L_edge_mid * sin(phi_cross));
fprintf(' %-30s : %8.4f in\n', 'Vertical Drop (Bike COG)', delta_bike_rack);
fprintf('%s\n', repmat('=',1,50));


%% 5. DYNAMIC ANALYSIS ("Pothole Test")

G_factor = 4.0; 
P_dynamic = P_load * G_factor;

fprintf('\n%s\n', repmat('=',1,50));
fprintf(' DYNAMIC ANALYSIS (%.1f G Pothole Load)\n', G_factor);
fprintf('%s\n', repmat('-',1,50));
fprintf(' %-30s : %8.2f lbs\n', 'Dynamic Load', P_dynamic);
fprintf('\n %-20s | %-12s | %-8s \n', 'COMPONENT', 'F.O.S.', 'STATUS');
fprintf('%s\n', repmat('-',1,50));

% Check 1: Weld (Dynamic)
fos_weld_dyn = Fmax_weld / P_dynamic;
if P_dynamic > Fmax_weld
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-8s\n', 'Weld Shear', fos_weld_dyn, status);

% Check 2: Crossbar Torsion (Dynamic)
shear_stress_dyn = shear_stress_cross * G_factor;
fos_bar_dyn = allowable_shear / shear_stress_dyn;
if shear_stress_dyn > allowable_shear
    status = '>> FAIL <<';
else
    status = 'PASS';
end
fprintf(' %-20s | %12.2f | %-8s\n', 'Crossbar Torsion', fos_bar_dyn, status);

% Check 3: Receiver Bending Stress (Dynamic)
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

% Deflection Check (Dynamic)
delta_dynamic_total_edge = delta_total_edge * G_factor;
delta_dynamic_drop = delta_bike_rack * G_factor;

fprintf(' %-30s : %8.4f in\n', 'Dynamic Tip Drop', delta_dynamic_total_edge);
fprintf('      > Due to Bend    : %8.4f in\n', delta_rec * G_factor);
fprintf('      > Due to Twist   : %8.4f in\n', (L_edge_mid * sin(phi_cross)) * G_factor);
fprintf(' %-30s : %8.4f in\n', 'Dynamic COG Drop', delta_dynamic_drop);
fprintf('%s\n', repmat('=',1,50));