close all
clear all
clc

set(0, 'DefaultFigureWindowStyle', 'docked')

% Input Data

D = 0.1;
Wb = 0.75;
rho = 998.23;
nu = 1.0e-6;
L = 3.0;
mu = nu * rho;
R = D / 2;

% Calculate Reynold's number

Re = D * Wb / nu;
fprintf('Reynolds Number = %.2f\n', Re)

% Claculate entrance length based on 2 formulas

Le1 = 40 * D;
Le2 = 60 * D;

fprintf('Entrance Length based on first formula = %.2f\n', Le1)
fprintf('Entrance Length based on second formula = %.2f\n', Le2)

% import mesh and extract mesh data

mesh = load('main.mat');
YPLS1 = mesh.YPLS;
YPLS1 = squeeze(YPLS1);
YPLS1 = YPLS1(31, :);
Z = mesh.Z_C;
Y = (mesh.Y_C)';
NZ = mesh.NZ;
NY = mesh.NY;
W = squeeze(mesh.W1);
P = squeeze(mesh.P1);
DWDY = squeeze(mesh.DWDY);
DPDZ = (P(1, 500) - P(1, 600)) / (Z(500) - Z(600));
ENUT = squeeze(mesh.ENUT);
mu_turb = rho * ENUT;
tau_re = - mu_turb .* DWDY;
tau_visc = - mu * DWDY;
tau_tot = tau_visc + tau_re;
tau_wall = rho * squeeze(mesh.STRS);
tau_wall = tau_wall(31, :);
W_tau = sqrt(tau_wall / rho);
WPLS = W ./ W_tau;
y_dummy = Y.* ones(1, 700);
YPLS = W_tau .* y_dummy / nu;
K = squeeze(mesh.KE);
EPS = squeeze(mesh.EP);

%%
% Checking 30 < YPLS < 130

figure('Color','w'); hold on; box on
plot(Z, YPLS1, 'LineWidth', 4)
xlabel('Axial Distance (m)', 'FontSize', 12, 'FontWeight','bold')
ylabel('Y+', 'FontSize', 12, 'FontWeight','bold')
ylim([45 75])
grid on

%  Velocity and tau_wall at 2.0m, 3.0m, 4.0m, 5.0m, 6.0m to verify fully developed flow

[~, i] = min(abs(Z - 2.0));
W_2000 = W(:, i);
tau_wall_2000 = tau_wall(i);

[~, i] = min(abs(Z - 3.0));
W_3000 = W(:, i);
tau_wall_3000 = tau_wall(i);

[~, i] = min(abs(Z - 4.0));
W_4000 = W(:, i);
tau_wall_4000 = tau_wall(i);

[~, i] = min(abs(Z - 5.0));
W_5000 = W(:, i);
tau_wall_5000 = tau_wall(i);

[~, i] = min(abs(Z - 6.0));
W_6000 = W(:, i);
tau_wall_6000 = tau_wall(i);


figure('Color','w'); hold on; box on
colors = lines(5);

plot(W_2000, Y, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(W_3000, Y, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(W_4000, Y, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(W_5000, Y, '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(W_6000, Y, '--','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))

xlabel('Velocity W (m/s)','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius R (m)','FontSize',12,'FontWeight','bold')
legend({'2.0 m','3.0 m','4.0 m','5.0 m', '6.0 m'}, 'Location','best', 'FontSize',11)

figure('Color','w'); hold on; box on
colors = lines(5);

plot(2.0,  tau_wall_2000,  'd','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(1,:),'MarkerFaceColor',colors(1,:))
plot(3.0,  tau_wall_3000, 'o','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(2,:),'MarkerFaceColor',colors(2,:))
plot(4.0, tau_wall_4000, 's','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(3,:),'MarkerFaceColor',colors(3,:))
plot(5.0, tau_wall_5000, '^','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(4,:),'MarkerFaceColor',colors(4,:))
plot(6.0,  tau_wall_6000, 'v','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(5,:),'MarkerFaceColor',colors(5,:))

xlim([1.0 7.0])
ylim([1 1.5])
grid on

xlabel('Axial Distance (m)','FontSize',12,'FontWeight','bold')
ylabel('Wall Shear Stress (Pa)','FontSize',12,'FontWeight','bold')
legend({'2.0 m','3.0 m','4.0 m','5.0 m', '6.0 m'}, 'Location','best', 'FontSize',11)

%%
% Pressure Distribution

Pressure = P(1,:);

% Friction factor

f_cfd = 4 * tau_wall / (0.5 * rho * Wb^2);
f_cfd_6000 = f_cfd(600);
f_blazius = 0.316 * Re^(-0.25);
f_blazius_plot = f_blazius * ones(1, 700);

fprintf('Friction factor value based on simulation = %.4f\n', f_cfd_6000)
fprintf('Friction factor value based on Blazius formula = %.4f\n', f_blazius)

figure('Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

% --- Pressure plot ---
nexttile
plot(Z, Pressure, 'LineWidth', 3)
box on
ylabel('Pressure (Pa)','FontSize',12,'FontWeight','bold')

% --- Friction factor plot ---
nexttile
plot(Z, f_cfd, 'LineWidth', 3)
box on
xlabel('Axial Distance (m)','FontSize',12,'FontWeight','bold')
ylabel('Friction Factor f','FontSize',12,'FontWeight','bold')
ylim([0.014 0.034])

%%
% Velocity profile in dimensionless form

[~, i] = min(abs(Z - 6.0));

WPLS_6000 = WPLS(:, i); 
YPLS_6000 = YPLS(:, i);
YPLS_6000 = flip(YPLS_6000);

figure('Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
loglog(WPLS_6000, YPLS_6000, 'LineWidth', 4)
xlabel('Velocity W+','FontSize',12,'FontWeight','bold')
ylabel('Wall Distance Y+','FontSize',12,'FontWeight','bold')

nexttile
plot(W_6000/Wb, Y/R, 'LineWidth', 4)
xlabel('Velocity W / Wb','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')

%%
% Total shear stress and its components

[~, i] = min(abs(Z - 6.0));
tau_visc_6000 = tau_visc(:, i);
tau_re_6000 = tau_re(:,i);
tau_tot_6000 = tau_tot(:, i);

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
plot(Y/R, tau_visc_6000 / tau_wall_6000, 'LineWidth', 3)
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Viscous Shear Stress / Tau w','FontSize',12,'FontWeight','bold')

nexttile
plot(Y/R, tau_re_6000 / tau_wall_6000, 'LineWidth', 3)
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Reynolds  Stress / Tau w','FontSize',12,'FontWeight','bold')

nexttile
hold on
plot(Y/R, tau_tot_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y/R, tau_visc_6000 / tau_wall_6000, 'LineStyle', '-','LineWidth', 3)
plot(Y/R, tau_re_6000 / tau_wall_6000, 'LineStyle', '--','LineWidth', 3)
ylim([0 1])
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Shear Stress Tau / Tau w','FontSize',12,'FontWeight','bold')
legend({'Tau Tot', 'Tau Visc', 'Tau Re'}, 'Location','best', 'FontSize',11)

%%
% Turbulent Kinetic Energy

[~, i] = min(abs(Z - 6.0));
K_6000 = K(:, i);

% Turbulent Dissipation Rate

[~, i] = min(abs(Z - 6.0));
EPS_6000 = EPS(:, i);

figure('Color','w')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
plot(Y/R, K_6000 / (Wb^2), 'LineWidth', 4)
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Turbulent Kinetic Energy K / Wb^2','FontSize',12,'FontWeight','bold')

nexttile
plot(Y/R, EPS_6000 * R / (Wb^3), 'LineWidth', 4)
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Turbulent Dissipation Rate EPS R / Wb^3','FontSize',12,'FontWeight','bold')

%%
% Grid Independence study

mesh1 = load('mesh1.mat');
mesh2 = load('mesh2.mat');
mesh3 = load('mesh3.mat');

W1 = mesh1.W1;
W1 = squeeze(W1);
P1 = mesh1.P1;
P1 = squeeze(P1);
NY1 = mesh1.NY;
NZ1 = mesh1.NZ;
tau_wall1 = rho * squeeze(mesh1.STRS);
tau_wall1 = tau_wall1(NY1, :);
Y1 = mesh1.Y_C;
Z1 = mesh1.Z_C;

W2 = mesh2.W1;
W2 = squeeze(W2);
P2 = mesh2.P1;
P2 = squeeze(P2);
NY2 = mesh2.NY;
NZ2 = mesh2.NZ;
tau_wall2 = rho * squeeze(mesh2.STRS);
tau_wall2 = tau_wall2(NY2, :);
Y2 = mesh2.Y_C;
Z2 = mesh2.Z_C;

W3 = mesh3.W1;
W3 = squeeze(W3);
P3 = mesh3.P1;
P3 = squeeze(P3);
NY3 = mesh3.NY;
NZ3 = mesh3.NZ;
tau_wall3 = rho * squeeze(mesh3.STRS);
tau_wall3 = tau_wall3(NY3, :);
Y3 = mesh3.Y_C;
Z3 = mesh3.Z_C; 

% Extracting values at 6.0 m from inlet and at centerline

[~, i] = min(abs(Z1 - 6.0));
W1_6000 = W1(:, i);
W1_loc = W1(1, i);
P1_loc = P1(1, i);
tau_wall1_6000 = tau_wall1(i);

[~, i] = min(abs(Z2 - 6.0));
W2_6000 = W2(:, i);
W2_loc = W2(1, i);
P2_loc = P2(1, i);
tau_wall2_6000 = tau_wall2(i);

[~, i] = min(abs(Z3 - 6.0));
W3_6000 = W3(:, i);
W3_loc = W3(1, i);
P3_loc = P3(1, i);
tau_wall3_6000 = tau_wall3(i);

[~, i] = min(abs(Z - 6.0));
W_6000 = W(:, i);
W_loc = W(1, i);
P_loc = P(1, i);
tau_wall_6000 = tau_wall(i);

figure('Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact')

nexttile
hold on; box on
colors = lines(4);

plot(NZ1, P1_loc, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ2, P2_loc, 'x', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', colors(2,:))

plot(NZ, P_loc, 's', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ3, P3_loc, 'd', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

xlim([350 850])
ylim([40 70])
xlabel('Number of Elements', 'FontSize', 12, 'FontWeight','bold')
ylabel('Pressure (Pa)', 'FontSize', 12, 'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)
grid on

nexttile
hold on; box on

plot(NZ1, W1_loc, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ2, W2_loc, 'x', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', colors(2,:))

plot(NZ, W_loc, 's', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ3, W3_loc, 'd', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

xlim([350 850])
ylim([0.7 1.2])
xlabel('Number of Elements', 'FontSize', 12, 'FontWeight','bold')
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)
grid on

nexttile
hold on; box on

plot(NZ1, tau_wall1_6000, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ2, tau_wall2_6000, 'x', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', colors(2,:))

plot(NZ, tau_wall_6000, 's', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ3, tau_wall3_6000, 'd', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

xlim([350 850])
ylim([1.2 1.5])
xlabel('Number of Elements', 'FontSize', 12, 'FontWeight','bold')
ylabel('Wall Shear Stress (Pa)', 'FontSize', 12, 'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)
grid on

figure('Color','w'); hold on; box on
colors = lines(4);

plot(W1_6000, Y1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(W2_6000, Y2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(W_6000, Y, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(W3_6000, Y3, '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

xlabel('Velocity W (m/s)','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius R (m)','FontSize',12,'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)

%%
% Validation

val_data = load('exp_data_Torgbersen.mat');
velocity_series1 = val_data.W_profile_series1;
w_series1 = velocity_series1(:,2);
y_series1 = velocity_series1(:,1);

velocity_series2 = val_data.W_profile_series2;
w_series2 = velocity_series2(:,2);
y_series2 = velocity_series2(:,1);

KE_series = val_data.k_profile;
k_val = KE_series(:, 2);
y_k = KE_series(:, 1);

epsilon_series1 = val_data.eps_profile_series1;
eps_series1 = epsilon_series1(:, 2);
y_eps_1 = epsilon_series1(:, 1);

epsilon_series2 = val_data.eps_profile_series2;
eps_series2 = epsilon_series2(:, 2);
y_eps_2 = epsilon_series2(:, 1);

epsilon_series3 = val_data.eps_profile_series3;
eps_series3 = epsilon_series3(:, 2);
y_eps_3 = epsilon_series3(:, 1);

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
hold on
plot(y_series1, w_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_series2, w_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(Y/R, W_6000/Wb, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))

ylabel('Velocity W / Wb','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Series 1','Series 2','CFD'}, 'Location','best', 'FontSize',11)


nexttile
hold on
plot(y_k, k_val, 'LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(Y/R, K_6000 / (Wb^2), '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))

ylabel('Turb. Kinetic Energy K / Wb^2','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Series 1','CFD'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(y_eps_1, eps_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_eps_2, eps_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(y_eps_3, eps_series3, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y/R, EPS_6000 * R / (Wb^3), '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

ylabel('Dissip. of Kinetic Energy EPS R / Wb^3','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Series 1','Series 2', 'Series 3','CFD'}, 'Location','best', 'FontSize',11)

%%
% Friction Factor

function f = prandt(Reb)
    prandt_eq = @(f) 1/sqrt(f) + 2 * log10(2.51/(Reb * sqrt(f)));
    f0 = 0.03;
    f = fzero(prandt_eq, f0);
end

function f = Haaland(Reb)
    Haaland_eq = @(f) 1/sqrt(f) + 1.8 * log10(6.9/(Reb));
    f0 = 0.03;
    f = fzero(Haaland_eq, f0);
end

f_prandt = prandt(Re);
f_haaland = Haaland(Re);
f_moody = 0.019;
fprintf('Friction Factor obtained by Prandt equation = %.4f\n', f_prandt)
fprintf('Friction Factor obtained by Haaland equation = %.4f\n', f_haaland)
fprintf('Friction Factor obtained by Moody chart = %.4f\n', f_moody)

%%
% Velocity using Nikuradse

n_prandt = f_prandt^(-0.5);
n_haaland = f_haaland^(-0.5);
n_moody = 0.019^(-0.5);

w_prandt = (1 / (2 * n_prandt^2)) * (2 * n_prandt + 1) * (n_prandt + 1) * (1 - Y / R).^(1 / n_prandt);
w_haaland = (1 / (2 * n_haaland^2)) * (2 * n_haaland + 1) * (n_haaland + 1) * (1 - Y / R).^(1 / n_haaland);
w_moody = (1 / (2 * n_moody^2)) * (2 * n_moody + 1) * (n_moody + 1) * (1 - Y / R).^(1 / n_moody);

figure('Color','w'); hold on
colors = lines(6);

plot(Y/R, w_prandt, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(Y/R, w_haaland, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(Y/R, w_moody, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(y_series1, w_series1, '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(y_series2, w_series2, '-k','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))
plot(Y/R, W_6000 / Wb, '--','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(6,:))

ylabel('Velocity W / Wb','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Prandt', 'Haaland', 'Moody', 'Series 1','Series 2', 'CFD'}, 'Location','best', 'FontSize',11)

%%
% Velocity using Launder and Spalding

k = 0.41;
E = 8.6;
indices = find(YPLS_6000 <= 200);
YPLS_filtered = YPLS_6000(indices);
WPLS_filtered = WPLS_6000(indices);

% calculate_w_plus = @(y) (y < 11.6) .* y + (y > 11.6 & y <= 130) .* (1 / k * log(E * y));
calculate_w_plus = @(y) (1 / k * log(E * y));

WPLS_LS = calculate_w_plus(YPLS_filtered);

figure('Color','w'); hold on
colors = lines(2);

loglog(WPLS_filtered, YPLS_filtered, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
loglog(WPLS_LS, YPLS_filtered, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))

xlabel('Velocity W+','FontSize',12,'FontWeight','bold')
ylabel('Wall Distance Y+','FontSize',12,'FontWeight','bold')
legend({'CFD', 'Launder and Spalding'}, 'Location','best', 'FontSize',11)

%%
% Turbulent Kinetic Energy

k_prandt = (f_prandt / 8) * (1 + (2/3) * (Y / R) + (10/3) * (Y / R).^3);
k_haaland = (f_haaland / 8) * (1 + (2/3) * (Y / R) + (10/3) * (Y / R).^3);
k_moody = (f_moody / 8) * (1 + (2/3) * (Y / R) + (10/3) * (Y / R).^3);

figure('Color','w'); hold on
colors = lines(5);

plot(Y/R, K_6000 / (Wb^2), '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_k, k_val, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(Y/R, k_moody, '-k','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))
plot(Y/R, k_prandt, '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y/R, k_haaland, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))

xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Turb. Kinetic Energy K / Wb^2','FontSize',12,'FontWeight','bold')
legend({'CFD', 'Torbergsen', 'Moody', 'Prandt', 'Haaland'}, 'Location','best', 'FontSize',11)

%%
% Dissipation Rate

lm = R * (0.14 - 0.08 * (Y / R).^2 - 0.06 * (Y / R).^4);

eps_prandt = 0.1643 * (k_prandt).^(3/2) ./ lm;
eps_haaland = 0.1643 * (k_haaland).^(3/2) ./ lm;
eps_moody = 0.1643 * (k_moody).^(3/2) ./ lm;


figure('Color','w'); hold on
colors = lines(7);

plot(Y/R, EPS_6000 * R / (Wb^3), '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(y_eps_1, eps_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_eps_2, eps_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(y_eps_3, eps_series3, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y/R, eps_moody, '-','LineWidth', 2.5, 'Color', colors(5,:), 'Marker','none')
plot(Y/R, eps_prandt, '--o','LineWidth', 2.5, 'Color', colors(6,:), 'Marker','none')
plot(Y/R, eps_haaland, '--d','LineWidth', 2.5, 'Color', colors(7,:), 'Marker','none')

ylim([0 0.1])

ylabel('Dissip. of Kinetic Energy EPS R / Wb^3','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'CFD', 'Series 1','Series 2', 'Series 3', 'Moody', 'Prandt', 'Haaland'}, 'Location','best', 'FontSize',11)

%%
% Influence of Turbulence Model

% Import Data of RNG

RNG = load('RNG.mat');


Z_RNG = RNG.Z_C;
Y_RNG = (RNG.Y_C)';
NZ_RNG = RNG.NZ;
NY_RNG = RNG.NY;
W_RNG = squeeze(RNG.W1);
P_RNG = squeeze(RNG.P1);
DWDY_RNG = squeeze(RNG.DWDY);
ENUT_RNG = squeeze(RNG.ENUT);
mu_turb_RNG = rho * ENUT_RNG;
tau_re_RNG = - mu_turb_RNG .* DWDY_RNG;
tau_visc_RNG = - mu * DWDY_RNG;
tau_tot_RNG = tau_visc_RNG + tau_re_RNG;
tau_wall_RNG = rho * squeeze(RNG.STRS);
tau_wall_RNG = tau_wall_RNG(31, :);
W_tau_RNG = sqrt(tau_wall_RNG / rho);
WPLS_RNG = W_RNG ./ W_tau_RNG;
y_dummy_RNG = Y_RNG.* ones(1, 700);
YPLS_RNG = W_tau_RNG .* y_dummy_RNG / nu;
K_RNG = squeeze(RNG.KE);
EPS_RNG = squeeze(RNG.EP);

[~, i] = min(abs(Z_RNG - 6.0));
W_RNG_6000 = W_RNG(:, i);
K_RNG_6000 = K_RNG(:, i);
EPS_RNG_6000 = EPS_RNG(:, i);
tau_wall_RNG_6000 = tau_wall_RNG(i);
tau_re_RNG_6000 = tau_re_RNG(:, i);
tau_visc_RNG_6000 = tau_visc_RNG(:, i);
tau_tot_RNG_6000 = tau_re_RNG_6000 + tau_visc_RNG_6000;

%%
% Import Data of Realisable

Real = load('Realisable.mat');

Z_Real = Real.Z_C;
Y_Real = (Real.Y_C)';
NZ_Real = Real.NZ;
NY_Real = Real.NY;
W_Real = squeeze(Real.W1);
P_Real = squeeze(Real.P1);
DWDY_Real = squeeze(Real.DWDY);
ENUT_Real = squeeze(Real.ENUT);
mu_turb_Real = rho * ENUT_Real;
tau_re_Real = - mu_turb_Real .* DWDY_Real;
tau_visc_Real = - mu * DWDY_Real;
tau_tot_Real = tau_visc_Real + tau_re_Real;
tau_wall_Real = rho * squeeze(Real.STRS);
tau_wall_Real = tau_wall_Real(31, :);
W_tau_Real = sqrt(tau_wall_Real / rho);
WPLS_Real = W_Real ./ W_tau_Real;
y_dummy_Real = Y_Real.* ones(1, 700);
YPLS_Real = W_tau_Real .* y_dummy_Real / nu;
K_Real = squeeze(Real.KE);
EPS_Real = squeeze(Real.EP);

[~, i] = min(abs(Z_Real - 6.0));
W_Real_6000 = W_Real(:, i);
K_Real_6000 = K_Real(:, i);
EPS_Real_6000 = EPS_Real(:, i);
tau_wall_Real_6000 = tau_wall_Real(i);
tau_re_Real_6000 = tau_re_Real(:, i);
tau_visc_Real_6000 = tau_visc_Real(:, i);
tau_tot_Real_6000 = tau_re_Real_6000 + tau_visc_Real_6000;

%%
% Validation

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
colors = lines(7);

nexttile
hold on
plot(y_series1, w_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_series2, w_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(Y/R, W_6000/Wb, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_RNG/R, W_RNG_6000/Wb, '-s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y_Real/R, W_Real_6000/Wb, '-s','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))

ylabel('Velocity W / Wb','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Series 1','Series 2','Standard K-e', 'RNG K-e', 'Realisable'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(y_k, k_val, 'LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(Y/R, K_6000 / (Wb^2), '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(Y_RNG/R, K_RNG_6000/ (Wb^2), '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, K_Real_6000/ (Wb^2), '-s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

ylabel('Turb. Kinetic Energy K / Wb^2','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Series 1','Standard K-e', 'RNG K-e', 'Realisable'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(y_eps_1, eps_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_eps_2, eps_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(y_eps_3, eps_series3, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y/R, EPS_6000 * R / (Wb^3), '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y_RNG/R, EPS_RNG_6000 * R / (Wb^3), '-s','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))
plot(Y_Real/R, EPS_Real_6000 * R / (Wb^3), '-s','LineWidth', 2.5, 'Color', colors(6,:), 'MarkerSize',5, 'MarkerFaceColor', colors(6,:))

ylim([0 0.002])

ylabel('Dissip. of Kinetic Energy EPS R / Wb^3','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
legend({'Series 1','Series 2', 'Series 3','Standard K-e', 'RNG K-e', 'Realisable'}, 'Location','best', 'FontSize',11)

%%
% Shear Stress

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
colors = lines(6);

nexttile
hold on
plot(Y/R, tau_visc_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y_RNG/R, tau_visc_RNG_6000 / tau_wall_RNG_6000, '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, tau_visc_Real_6000 / tau_wall_Real_6000, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Viscous Shear Stress / Tau w','FontSize',12,'FontWeight','bold')
legend({'Standard K-e', 'RNG K-e', 'Realisable'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(Y/R, tau_re_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y_RNG/R, tau_re_RNG_6000 / tau_wall_RNG_6000, '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, tau_re_Real_6000 / tau_wall_Real_6000, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Reynolds  Stress / Tau w','FontSize',12,'FontWeight','bold')
legend({'Standard K-e', 'RNG K-e', 'Realisable'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(Y/R, tau_tot_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y_RNG/R, tau_tot_RNG_6000 / tau_wall_RNG_6000, '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, tau_tot_Real_6000 / tau_wall_Real_6000, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

ylim([0 1])
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Shear Stress Tau / Tau w','FontSize',12,'FontWeight','bold')
legend({'Standard K-e', 'RNG K-e', 'Realisable'}, 'Location','best', 'FontSize',11)

%%
% Low Reynold's

LRE = load('LowRe.mat');
YPLS1_LRE = LRE.YPLS;
YPLS1_LRE = squeeze(YPLS1_LRE);
NY_LRE = LRE.NY;
YPLS1_LRE = YPLS1_LRE(NY_LRE, :);

Z_LRE = LRE.Z_C;
Y_LRE = (LRE.Y_C)';
NZ_LRE = LRE.NZ;
W_LRE = squeeze(LRE.W1);
P_LRE = squeeze(LRE.P1);
DWDY_LRE = squeeze(LRE.DWDY);
ENUT_LRE = squeeze(LRE.ENUT);
mu_turb_LRE = rho * ENUT_LRE;
tau_re_LRE = - mu_turb_LRE .* DWDY_LRE;
tau_visc_LRE = - mu * DWDY_LRE;
tau_tot_LRE = tau_visc_LRE + tau_re_LRE;
% tau_wall_LRE = rho * squeeze(LRE.STRS);
tau_wall_LRE = - mu * DWDY_LRE;
tau_wall_LRE = tau_wall_LRE(NY_LRE, :);
W_tau_LRE = sqrt(tau_wall_LRE / rho);
WPLS_LRE = W_LRE ./ W_tau_LRE;
y_dummy_LRE = Y_LRE.* ones(1, NZ_LRE);
YPLS_LRE = W_tau_LRE .* y_dummy_LRE / nu;
K_LRE = squeeze(LRE.KE);
EPS_LRE = squeeze(LRE.EP);

[~, i] = min(abs(Z_LRE - 6.0));
W_LRE_6000 = W_LRE(:, i);
K_LRE_6000 = K_LRE(:, i);
EPS_LRE_6000 = EPS_LRE(:, i);
tau_wall_LRE_6000 = tau_wall_LRE(i);
tau_re_LRE_6000 = tau_re_LRE(:, i);
tau_visc_LRE_6000 = tau_visc_LRE(:, i);
tau_tot_LRE_6000 = tau_re_LRE_6000 + tau_visc_LRE_6000;

%%
% Checking 3YPLS ~ 1

figure('Color','w'); hold on; box on
plot(Z_LRE, YPLS1_LRE, 'LineWidth', 4)
xlabel('Axial Distance (m)', 'FontSize', 12, 'FontWeight','bold')
ylabel('Y+', 'FontSize', 12, 'FontWeight','bold')
%ylim([45 75])
grid on

%%
% Validation

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
colors = lines(7);

nexttile
hold on
plot(y_series1, w_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_series2, w_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(Y/R, W_6000/Wb, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
% plot(Y_RNG/R, W_RNG_6000/Wb, '-s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
% plot(Y_Real/R, W_Real_6000/Wb, '-s','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))
plot(Y_LRE/R, W_LRE_6000/Wb, '-s','LineWidth', 2.5, 'Color', colors(6,:), 'MarkerSize',5, 'MarkerFaceColor', colors(6,:))

ylabel('Velocity W / Wb','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
% legend({'Series 1','Series 2','Standard K-e', 'RNG K-e', 'Realisable', 'Low Reynolds'}, 'Location','best', 'FontSize',11)
legend({'Series 1','Series 2', 'Standard K-e', 'Low Reynolds'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(y_k, k_val, 'LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(Y/R, K_6000 / (Wb^2), '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
% plot(Y_RNG/R, K_RNG_6000/ (Wb^2), '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
% plot(Y_Real/R, K_Real_6000/ (Wb^2), '-s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y_LRE/R, K_LRE_6000/ (Wb^2), '-s','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))

ylabel('Turb. Kinetic Energy K / Wb^2','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
% legend({'Series 1','Standard K-e', 'RNG K-e', 'Realisable', 'Low Reynolds'}, 'Location','best', 'FontSize',11)
legend({'Series 1', 'Standard K-e', 'Low Reynolds'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(y_eps_1, eps_series1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(y_eps_2, eps_series2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(y_eps_3, eps_series3, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y/R, EPS_6000 * R / (Wb^3), '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
% plot(Y_RNG/R, EPS_RNG_6000 * R / (Wb^3), '-s','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))
% plot(Y_Real/R, EPS_Real_6000 * R / (Wb^3), '-s','LineWidth', 2.5, 'Color', colors(6,:), 'MarkerSize',5, 'MarkerFaceColor', colors(6,:))
plot(Y_LRE/R, EPS_LRE_6000 * R / (Wb^3), '-s','LineWidth', 2.5, 'Color', colors(7,:), 'MarkerSize',5, 'MarkerFaceColor', colors(7,:))

%ylim([0 0.002])

ylabel('Dissip. of Kinetic Energy EPS R / Wb^3','FontSize',12,'FontWeight','bold')
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
% legend({'Series 1','Series 2', 'Series 3','Standard K-e', 'RNG K-e', 'Realisable', 'Low Reynolds'}, 'Location','best', 'FontSize',11)
legend({'Series 1','Series 2', 'Series 3', 'Standard K-e', 'Low Reynolds'}, 'Location','best', 'FontSize',11)

%%
% Shear Stress

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
colors = lines(6);

nexttile
hold on
plot(Y/R, tau_visc_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y_RNG/R, tau_visc_RNG_6000 / tau_wall_RNG_6000, '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, tau_visc_Real_6000 / tau_wall_Real_6000, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y_LRE/R, tau_visc_LRE_6000 / tau_wall_LRE_6000, '--o','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))

xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Viscous Shear Stress / Tau w','FontSize',12,'FontWeight','bold')
legend({'Standard K-e', 'RNG K-e', 'Realisable', 'Low Reynolds'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(Y/R, tau_re_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y_RNG/R, tau_re_RNG_6000 / tau_wall_RNG_6000, '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, tau_re_Real_6000 / tau_wall_Real_6000, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y_LRE/R, tau_re_LRE_6000 / tau_wall_LRE_6000, '--o','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))


xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Reynolds  Stress / Tau w','FontSize',12,'FontWeight','bold')
legend({'Standard K-e', 'RNG K-e', 'Realisable', 'Low Reynolds'}, 'Location','best', 'FontSize',11)

nexttile
hold on
plot(Y/R, tau_tot_6000 / tau_wall_6000, 'LineWidth', 3)
plot(Y_RNG/R, tau_tot_RNG_6000 / tau_wall_RNG_6000, '-s','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(Y_Real/R, tau_tot_Real_6000 / tau_wall_Real_6000, '--o','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(Y_LRE/R, tau_tot_LRE_6000 / tau_wall_LRE_6000, '--o','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))

%ylim([0 1])
xlabel('Pipe Radius r / R','FontSize',12,'FontWeight','bold')
ylabel('Shear Stress Tau / Tau w','FontSize',12,'FontWeight','bold')
legend({'Standard K-e', 'RNG K-e', 'Realisable', 'Low Reynolds'}, 'Location','best', 'FontSize',11)