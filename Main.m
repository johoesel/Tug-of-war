%% Antigen extraction under Force
% This script runs the functions antigen_extract_dynamics.m, dependent on
% an energy landscape described by Potentials.m to understand the escape
% probability of an antigen molecule, given a force applied to the BCR-bond
%--------15/04/24--------
%Josephine Hoesel
clear;clc;close all
%%
%Physical Constants
gamma_a = 15e-4;
gamma_b = 15e-4;
gamma_r = 15e-4;

KB = 1.38e-23;
T = 273;

%Energy landscape parameters
x_aF = 3e-9; %APC-antigen rupture length
x_bF = 2e-9; %BCR-antigen rupture length
DeltaG_a = 10*KB*T; %Energy barrier
DeltaG_b = DeltaG_a;

%Initialise bond lengths (nm)
x_a = -2e-9:0.01e-9:4e-9;
x_b = x_a;
theta = 0;
r = 0;
F = 0;
x0 = 0;
y0 = 0;

[~,~,~,~,U_a,U_b] = Potentials(x_a,x_b,x0,y0,x_aF,x_bF,DeltaG_a,DeltaG_b,theta,r);

U = U_a + U_b - F * (x_a + x_b);

%Visualise the energy potential (1d)
figure(1)
plot(x_a*10^9,U/(KB*T),'k','LineWidth',1.5)
hold on

F = 20e-12;
U = U_a + U_b - F * (x_a + x_b);

plot(x_a*10^9,U/(KB*T),'--r','LineWidth',1.5)
yline(0)
xline(0)
xlabel('Bond extension (nm)','FontSize',14)
ylabel('Bond Free Energy (K_B T)','FontSize',14)

legend('Force = 0pN','Force = 40pN','FontSize',12)
exportgraphics(gcf,'Cubic_well.png')
%% Solving the SDE to find the antigen solution path and averaging 
%extraction probability over the number of runs, Runs > 1000

%Simulation Parameters
Tspan = 100;
dt = 0.001; %step-size
Runs = 1e4;
%Runs = 1;
%Testing the effect of force/rupture length
Forces = linspace(20e-12,200e-12,15);
D = 0;
%Forces = 80e-12;
x_aF = [0.5e-9,1e-9,1.5e-9,2e-9,3e-9];


Extraction_Probabilities = zeros(length(x_aF),length(Forces));
tic
for i = 1:length(x_aF)
    [Extraction_Probability,bond_coordinates,t_vec] = antigen_extract_dynamics(Tspan,dt,x_aF(i),x_bF,DeltaG_a,DeltaG_b,gamma_a,gamma_b,gamma_r,KB,T,Forces,D,Runs);
    Extraction_Probabilities(i,:) = Extraction_Probability;
end
toc
save("Extraction_probabilities_x_b_1p5_gabr_15e-4_D_0.mat","Extraction_Probabilities")
DF_Extract_Prob = array2table(Extraction_Probabilities,'VariableNames', 'F = '+ string(Forces*1e12) + 'pN', 'RowNames', string(x_aF*10^9)+'nm');
writetable(DF_Extract_Prob, (['Extraction_probabilities_x_b_1p5_gabr_15e-4_D_0', '.csv']),'WriteRowNames',true)

%% Visualising the brownian simulation on the 2d landscape
%Make Energy 2-dimensional
U_a_matrix = repmat(U_a,length(U_a),1);
U_b_matrix = repmat(U_b',1,length(U_b));
F = 0;
Force_matrix  = F * (repmat(x_a,length(x_a),1) + repmat(x_b',1,length(x_b)));
Total_Energy = U_a_matrix + U_b_matrix - Force_matrix;

%Plot one random trajectory
figure(2)

contourf(x_a(50:550)*10^9,x_b(50:550)*10^9,Total_Energy(50:550,50:550)/(KB*T),20, 'LineStyle','none')
shading interp
hold on
%set(gca,'ColorScale','log');
c = colorbar;
c.Label.String = 'Energy [K_BT]';
%plot(bond_coordinates(:,1)*10^9,bond_coordinates(:,2)*10^9,'Linewidth',1.5,color='k')
hold on
plot(0,0,'o','MarkerSize',8,MarkerFaceColor ='white')
xline(x_aF(end)*10^9, 'k--',{'x_a^F'})
yline(x_bF*10^9, 'k--',{'x_b^F'})
xlabel('APC-antigen bond extension, x_a (nm)')
ylabel('BCR-antigen bond extension, x_b (nm)')
exportgraphics(gcf,'random_trajectory.png')

%% Plot Force-dependent Extraction Probabilities
figure(3)
Forces = linspace(20,200,15);
x_aF = [0.5e-9,1e-9,1.5e-9,2e-9,3e-9];
for i = 1:length(x_aF)
    semilogy(Forces,Extraction_Probabilities(i,:), '-o','LineWidth',2)
    hold on
end

xlabel('Applied BCR Force (pN)')
ylabel('Extraction Probabilities, \eta')
legend("x_a^F = 0.5nm","x_a^F = 1.0nm","x_a^F = 1.5nm","x_a^F = 2.0nm","x_a^F = 3.0nm")
exportgraphics(gcf,'Extraction_vs_Force.png')




%% Testing the Effect of APC Mobility on the simulation

Tspan = 100;
dt = 0.001;
F = 0;
Runs = 1e4;
x_aF = x_bF;

U_a_matrix = repmat(U_a,length(U_a),1);
U_b_matrix = repmat(U_b',1,length(U_b));

%2.11 0.08

Forces = linspace(20e-12,200e-12,15);
Forces = 80e-12;
x_aF = [0.5e-9,1e-9,1.5e-9,2e-9,3e-9];
D = 0:0.2e-6:2.2e-6;


Extraction_Probabilities_mobility = zeros(length(x_aF),length(D));
final_Force_on_r = zeros(length(x_aF),length(D));
tic
for i = 1:length(x_aF)
    [Extraction_Probability_mobility,bond_coordinates,tvec,Theta,Fr] = antigen_extract_dynamics(Tspan,dt,x_aF(i),x_bF,DeltaG_a,DeltaG_b,gamma_a,gamma_b,gamma_r,KB,T,Forces,D,Runs);
    Extraction_Probabilities_mobility(i,:) = Extraction_Probability_mobility;
    final_Force_on_r(i,:) = Fr;
end

toc
save("Extraction_probabilities_mobility_x_b_1p5_gabr_15e-4_D_0_2.2em6_F80em12.mat","Extraction_Probabilities_mobility")
DF_Extract_Prob_m = array2table(Extraction_Probabilities_mobility,'VariableNames', 'D = '+ string(D*1e5), 'RowNames', string(x_aF*10^9)+'nm');
writetable(DF_Extract_Prob_m, (['Extraction_probabilities_mobility_x_b_1p5_gabr_15e-4_D_0_2.2em6_F80em12', '.csv']),'WriteRowNames',true)

Force_matrix  = F * (repmat(x_a,length(x_a),1) + repmat(x_b',1,length(x_b)));
Total_Energy = U_a_matrix + U_b_matrix - Force_matrix;

%Plot one random trajectory
figure(5)
contourf(x_a(50:550),x_b(50:550),Total_Energy(50:550,50:550)/(KB*T),20, 'LineStyle','none')
shading interp
hold on
%set(gca,'ColorScale','log');
c = colorbar;
c.Label.String = 'Energy [K_BT]';
plot(bond_coordinates(:,1),bond_coordinates(:,2),'Linewidth',1.5,color='k')
hold on
plot(0,0,'o','MarkerSize',8,MarkerFaceColor ='white')
xline(x_aF(end), 'k--',{'x_a^F'})
yline(x_bF, 'k--',{'x_b^F'})
xlabel('APC-antigen bond extension, x_a (nm)')
ylabel('BCR-antigen bond extension, x_b (nm)')
exportgraphics(gcf,'random_trajectory.png')


%% %% Plot Diffusion-dependent Extraction Probabilities
figure(6)
x_aF = [0.5e-9,1e-9,1.5e-9,2e-9,3e-9];
D = 0:0.2:2.2;
for i = 1:length(x_aF)
    semilogy(D,Extraction_Probabilities_mobility(i,:), '-o','LineWidth',2)
    hold on
end

xlabel('Diffusion Coefficient (\mu m^2 s^{-1})')
ylabel('Extraction Probabilities, \eta')
legend("x_a^F = 0.5nm","x_a^F = 1.0nm","x_a^F = 1.5nm","x_a^F = 2.0nm","x_a^F = 3.0nm")
exportgraphics(gcf,'Extraction_vs_Diffusion.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%