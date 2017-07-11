testOH = AcidBase();

pH = [1:0.1:14];
pK_OH = 16;
pK_OH2 = -2.4;

% Calculate scaled rate constant

label1 = { 'H_2O', 'OH^-', 'H_2PO_4^-', 'HPO_4^{2-}', 'PO_4^{3-}' };

[krate1_H2O] = testOH.rateReactionPerCatalyst(pK_OH, 'H2O', pH);
[krate1_OH] = testOH.rateReactionPerCatalyst(pK_OH, 'OH-', pH);
[krate1_H2PO4] = testOH.rateReactionPerCatalyst(pK_OH, 'H2PO4', pH);
[krate1_HPO4] = testOH.rateReactionPerCatalyst(pK_OH, 'HPO4', pH);
[krate1_PO4] = testOH.rateReactionPerCatalyst(pK_OH, 'PO4', pH);

label2 = { 'H_3O^+', 'H_2O', 'H_3PO_4', 'H_2PO_4^-', 'HPO_4^{2-}' };

[krate2_H3O] = testOH.rateReactionPerCatalyst2(pK_OH2, 'H3O+', pH);
[krate2_H2O] = testOH.rateReactionPerCatalyst2(pK_OH2, 'H2O', pH);
[krate2_H3PO4] = testOH.rateReactionPerCatalyst2(pK_OH2, 'H3PO4', pH);
[krate2_H2PO4] = testOH.rateReactionPerCatalyst2(pK_OH2, 'H2PO4', pH);
[krate2_HPO4] = testOH.rateReactionPerCatalyst2(pK_OH2, 'HPO4', pH);

label3 = { 'Water', 'PBS'};
[krate_Water] = testOH.rateReactionInWater2(pK_OH, pK_OH2, pH);
[krate_PBS] = testOH.rateReactionInPBS2(pK_OH, pK_OH2, pH);

% Set up for Bloch simulation

chemical_shift = 1.0 * 500;

T1_I = 4; % T1 of the abundant pool
T2_I = 0.06; % T2 of the abundant pool
T1_S = 4; % T1 of the solute pool
T2_S = 0.06; % T2 of the solute pool

testOH = testOH.ParametersForBlochMcConnell(chemical_shift, T1_I, T2_I, T1_S, T2_S)

% Anothre set up for Bloch simulation

concentration = 0.01; % in mol

w1 = 2.0 * pi * [100 0]; % the amplitudes of the saturating RF field in Hz
sat_time = 20; % the duration of the pre-saturation (sec)

freq_max = 5000;
freq_step = 10;

catalysts1 = {'H2O', 'OH-', 'H2PO4', 'HPO4', 'PO4'} ;
catalysts2 = {'H3O+', 'H2O', 'H3PO4', 'H2PO4', 'HPO4'} ;

tic

[MTRasym1_H2O, freq_offsets] = testOH.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'H2O', pK_OH, concentration);
[MTRasym1_OH, freq_offsets] = testOH.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'OH-', pK_OH, concentration);
[MTRasym1_H2PO4, freq_offsets] = testOH.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'H2PO4', pK_OH, concentration);
[MTRasym1_HPO4, freq_offsets] = testOH.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'HPO4', pK_OH, concentration);
[MTRasym1_PO4, freq_offsets] = testOH.MTRasym(freq_max, freq_step, w1, sat_time, pH, 'PO4', pK_OH, concentration);

[MTRasym2_H3O, freq_offsets] = testOH.MTRasym2(freq_max, freq_step, w1, sat_time, pH, 'H3O+', pK_OH, concentration, pK_OH2);
[MTRasym2_H2O, freq_offsets] = testOH.MTRasym2(freq_max, freq_step, w1, sat_time, pH, 'H2O', pK_OH, concentration, pK_OH2);
[MTRasym2_H3PO4, freq_offsets] = testOH.MTRasym2(freq_max, freq_step, w1, sat_time, pH, 'H3PO4', pK_OH, concentration, pK_OH2);
[MTRasym2_H2PO4, freq_offsets] = testOH.MTRasym2(freq_max, freq_step, w1, sat_time, pH, 'H2PO4', pK_OH, concentration, pK_OH2);
[MTRasym2_HPO4, freq_offsets] = testOH.MTRasym2(freq_max, freq_step, w1, sat_time, pH, 'HPO4', pK_OH, concentration, pK_OH2);

toc

tic
[MTRasym_Water, freq_offsets] = testOH.MTRasymWater2(freq_max, freq_step, w1, sat_time, pH, pK_OH, concentration, pK_OH2);
[MTRasym_PBS, freq_offsets] = testOH.MTRasymPBS2(freq_max, freq_step, w1, sat_time, pH, pK_OH, concentration, pK_OH2);

toc

%%
% Plots
ind = find(freq_offsets < testOH.chemical_shift, 1,'last') + 1;

figure('Position', [0 0 1500 800])
subplot(2,3,1)
semilogy(pH, krate1_H2O / 10^10, pH, krate1_OH / 10^10, ...
    pH, krate1_H2PO4 / 10^10, pH, krate1_HPO4 / 10^10, pH, krate1_PO4 / 10^10, ...
    'LineWidth', 3)

title(['Primary Alcohol, pK = 16'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label1, 'Location', 'northwest')
axis([1 14 10^-25 10^0])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,3,2)
semilogy(pH, krate2_H3O / 10^10, pH, krate2_H2O / 10^10, ...
    pH, krate2_H3PO4 / 10^10, pH, krate2_H2PO4 / 10^10, pH, krate2_HPO4 / 10^10, ...
    'LineWidth', 3)

title(['Primary Alcohol, pK = 16'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label2, 'Location','northeast')
axis([1 14 10^-25 10^0])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,3,3)
semilogy(pH, krate_Water / 10^10, pH, krate_PBS / 10^10, 'LineWidth', 3)

title(['Primary Alcohol, pK = 16'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label3, 'Location','north')
axis([1 14 10^-10 10^0])
set(gca,'FontSize',16,'XTick',1:14)


% MTRasym

subplot(2,3,4)
plot(pH, MTRasym1_H2O(ind,:), pH, MTRasym1_OH(ind,:), ...
    pH, MTRasym1_H2PO4(ind,:), pH, MTRasym1_HPO4(ind,:), pH, MTRasym1_PO4(ind,:), ...
    'LineWidth', 3)

title(['Primary Alcohol, pK = 16'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label1 ,'Location','west')
axis([1 14 0 0.02])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,3,5)
plot(pH, MTRasym2_H3O(ind,:), pH, MTRasym2_H2O(ind,:), ...
    pH, MTRasym2_H3PO4(ind,:), pH, MTRasym2_H2PO4(ind,:), pH, MTRasym2_HPO4(ind,:), ...
    'LineWidth', 3)

title(['Primary Alcohol, pK = 16'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label2 ,'Location','east')
axis([1 14 0 0.02])
set(gca,'FontSize',16,'XTick',1:14)

subplot(2,3,6)
plot(pH, MTRasym_Water(ind,:), pH, MTRasym_PBS(ind,:), 'LineWidth', 3)

title(['Primary Alcohol, pK = 16'],'FontSize',16)
xlabel('pH','FontSize',16)
ylabel('$\dot{\xi}$/($k_0$[HA])','FontSize',16,'Interpreter','Latex')
legend(label3 ,'Location','east')
axis([1 14 0 0.02])
set(gca,'FontSize',16,'XTick',1:14)
