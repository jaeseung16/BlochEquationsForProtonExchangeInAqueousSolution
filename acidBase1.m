testCr = AcidBase();
testGlu = AcidBase();

%
pH = [1:0.1:14];
pK_Cr = 11.02;
pK_Glu = 9.67;

%
% Calculate scaled rate constant

label1 = { 'H_2O', 'OH^-', 'H_2PO_4^-', 'HPO_4^{2-}', 'PO_4^{3-}' };

krate_H2O_Cr = testCr.rateReactionPerCatalyst(pK_Cr, 'H2O', pH);
krate_OH_Cr = testCr.rateReactionPerCatalyst(pK_Cr, 'OH-', pH);
krate_H2PO4_Cr = testCr.rateReactionPerCatalyst(pK_Cr, 'H2PO4', pH);
krate_HPO4_Cr = testCr.rateReactionPerCatalyst(pK_Cr, 'HPO4', pH);
krate_PO4_Cr = testCr.rateReactionPerCatalyst(pK_Cr, 'PO4', pH);

krate_H2O_Glu = testGlu.rateReactionPerCatalyst(pK_Glu, 'H2O', pH);
krate_OH_Glu = testGlu.rateReactionPerCatalyst(pK_Glu, 'OH-', pH);
krate_H2PO4_Glu = testGlu.rateReactionPerCatalyst(pK_Glu, 'H2PO4', pH);
krate_HPO4_Glu = testGlu.rateReactionPerCatalyst(pK_Glu, 'HPO4', pH);
krate_PO4_Glu = testGlu.rateReactionPerCatalyst(pK_Glu, 'PO4', pH);

label2 = { 'Protonated Amine in Water', 'Protonated Amine in PBS', ...
    'Guanidinium in Water', 'Guanidinium in PBS' };
krate_Water_Cr = testCr.rateReactionInWater(pK_Cr, pH);
krate_PBS_Cr = testCr.rateReactionInPBS(pK_Cr, pH);

krate_Water_Glu = testGlu.rateReactionInWater(pK_Glu, pH);
krate_PBS_Glu = testGlu.rateReactionInPBS(pK_Glu, pH);

%
% Plots of scaled rate constants

figure('Position',[0, 0, 1000, 1000])

subplot(2,2,1)
semilogy(pH, krate_H2O_Glu / 10^10, pH, krate_OH_Glu / 10^10, ...
    pH, krate_H2PO4_Glu / 10^10, pH, krate_HPO4_Glu / 10^10, pH, krate_PO4_Glu / 10^10, ...
    'LineWidth', 3)

title(['Protonated Amine, pK = 9.67'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label1 ,'Location','South')
axis([1 14 10^-25 10^0])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,2,2)
semilogy(pH, krate_H2O_Cr / 10^10, pH, krate_OH_Cr / 10^10, ...
    pH, krate_H2PO4_Cr / 10^10, pH, krate_HPO4_Cr / 10^10, pH, krate_PO4_Cr / 10^10, ...
    'LineWidth', 3)

title(['Guanidinium, pK = 11.02'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label1 ,'Location','South')
axis([1 14 10^-25 10^0])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,2,3)
semilogy(pH, krate_Water_Glu / 10^10, pH, krate_PBS_Glu / 10^10, ...
    pH, krate_Water_Cr / 10^10, pH, krate_PBS_Cr / 10^10, 'LineWidth', 3 )

title(['Summed scaled rates'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend('Protonated Amine in Water', 'Protonated Amine in PBS', ...
    'Guanidinium in Water', 'Guanidinium in PBS','Location','South')
axis([1 14 10^-25 10^0])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,2,4)
semilogy(pH, krate_Water_Glu / 10^10, pH, krate_PBS_Glu / 10^10, ...
    pH, krate_Water_Cr / 10^10, pH, krate_PBS_Cr / 10^10, 'LineWidth', 3 )

title(['Summed scaled rates (inlet)'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label2,'Location','South')
axis([1 7 10^-15 10^-5])
set(gca,'FontSize',16,'XTick',1:14)

%
% Set up for Bloch simulation

chemical_shift_Cr = 1.9 * 500;
chemical_shift_Glu = 3.0 * 500;

T1_I = 4; % T1 of the abundant pool
T2_I = 0.06; % T2 of the abundant pool
T1_S = 4; % T1 of the solute pool
T2_S = 0.06; % T2 of the solute pool

testCr = testCr.ParametersForBlochMcConnell(chemical_shift_Cr, T1_I, T2_I, T1_S, T2_S)
testGlu = testGlu.ParametersForBlochMcConnell(chemical_shift_Glu, T1_I, T2_I, T1_S, T2_S)

%
% Anothre set up for Bloch simulation
concentration = 0.01; % in mol
freq_offsets = -10000:10:10000; % the range of the frequency offsets in Hz
w1 = 2.0 * pi * [100 0]; % the amplitudes of the saturating RF field in Hz
sat_time = 20; % the duration of the pre-saturation (sec)

freq_max = 5000;
freq_step = 10;

%
catalysts = {'H2O', 'OH-', 'H2PO4', 'HPO4', 'PO4'} ;

tic

[MTRasym_H2O_Cr, freq_offsets] = testCr.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'H2O', pK_Cr, concentration);
[MTRasym_H2O_Glu, freq_offsets] = testGlu.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'H2O', pK_Glu, concentration);

[MTRasym_OH_Cr, freq_offsets] = testCr.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'OH-', pK_Cr, concentration);
[MTRasym_OH_Glu, freq_offsets] = testGlu.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'OH-', pK_Glu, concentration);

[MTRasym_H2PO4_Cr, freq_offsets] = testCr.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'H2PO4', pK_Cr, concentration);
[MTRasym_H2PO4_Glu, freq_offsets] = testGlu.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'H2PO4', pK_Glu, concentration);

[MTRasym_HPO4_Cr, freq_offsets] = testCr.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'HPO4', pK_Cr, concentration);
[MTRasym_HPO4_Glu, freq_offsets] = testGlu.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'HPO4', pK_Glu, concentration);

[MTRasym_PO4_Cr, freq_offsets] = testCr.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'PO4', pK_Cr, concentration);
[MTRasym_PO4_Glu, freq_offsets] = testGlu.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'PO4', pK_Glu, concentration);

[MTRasym_Water_Cr, freq_offsets] = testCr.MTRasymWater(freq_max, freq_step, w1, sat_time, pH, pK_Cr, concentration);
[MTRasym_Water_Glu, freq_offsets] = testGlu.MTRasymWater(freq_max, freq_step, w1, sat_time, pH, pK_Glu, concentration);

[MTRasym_PBS_Cr, freq_offsets] = testCr.MTRasymPBS(freq_max, freq_step, w1, sat_time, pH, pK_Cr, concentration);
[MTRasym_PBS_Glu, freq_offsets] = testGlu.MTRasymPBS(freq_max, freq_step, w1, sat_time, pH, pK_Glu, concentration);
toc

%

ind_Cr = find(freq_offsets < testCr.chemical_shift, 1,'last') + 1;
ind_Glu = find(freq_offsets < testGlu.chemical_shift, 1,'last') + 1;

figure('Position',[0 0 1500 500])
subplot(1,3,1)
plot(pH, MTRasym_H2O_Glu(ind_Glu,:), pH, MTRasym_OH_Glu(ind_Glu,:), ...
    pH, MTRasym_H2PO4_Glu(ind_Glu,:), pH, MTRasym_HPO4_Glu(ind_Glu,:), ...
    pH, MTRasym_PO4_Glu(ind_Glu,:), 'LineWidth', 3)

title(['Protonated Amine, pK = 9.67'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label1 ,'Location','East')
axis([1 14 0 0.15])
set(gca,'FontSize',16,'XTick',1:14)

subplot(1,3,2)
plot(pH, MTRasym_H2O_Cr(ind_Cr,:), pH, MTRasym_OH_Cr(ind_Cr,:), ...
    pH, MTRasym_H2PO4_Cr(ind_Cr,:), pH, MTRasym_HPO4_Cr(ind_Cr,:), ...
    pH, MTRasym_PO4_Cr(ind_Cr,:), 'LineWidth', 3)

title(['Guanidinium, pK = 11.02'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label1 ,'Location','East')
axis([1 14 0 0.15])
set(gca,'FontSize',16,'XTick',1:14)

subplot(1,3,3)
plot(pH, MTRasym_Water_Glu(ind_Glu,:), pH, MTRasym_PBS_Glu(ind_Glu,:), ...
    pH, MTRasym_Water_Cr(ind_Cr,:), pH, MTRasym_PBS_Cr(ind_Cr,:), 'LineWidth', 3)

title(['Combined reations'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label2,'Location','North')
axis([1 14 0 0.15])
set(gca,'FontSize',16,'XTick',1:14)
