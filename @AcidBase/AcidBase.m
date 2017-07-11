classdef AcidBase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pK_catalysis = containers.Map( {'H3O+','H2O','H3PO4','H2PO4','HPO4'}, ...
            [ -1.7, 15.7, 2.1, 7.2, 12.3] );
        
        pK_groups = containers.Map( ...
            {'OH2', 'primaryAlcohols','NH4+','NH3','amides','guanidinium','Glu','Cr','Amide'}, ...
            [-2.4, 16, 9.2, 38, 15, 13.6, 9.67, 11.02, 15] );
        
        k0 = 10^10; % maximum rate constant controlled by diffusion (mol-1 s-1)
        
        K1 = 10^(-2.1); % K of H3PO4
        K2 = 10^(-7.2); % K of H2PO4+
        K3 = 10^(-12.3); % K of HPO4^{2-}
        phosphate_concentration = 0.01; % Total concentration of inorganic phosphate (e.g. in PBS)
        
        chemical_shift % chemical shift in Hz
        T1_I % T1 of abundant pool in sec
        T2_I % T1 of abundant pool in sec
        T1_S % T1 of cest pool in sec
        T2_S % T2 of cest pool in sec
        
    end
    
    methods
        function obj = AcidBase()
            
        end
        
        % Methods related to chemical reactions
        
        [rateF] = rateConstant(obj, pK_donor, pK_acceptor)
        
        [phos1, phos2, phos3, phos4] = phosphateConcentration(obj, pH)
        
        % Reaction rate for reactions involving HA and A-
        function [krate] = rateReactionPerCatalyst(obj, pK_donor, catalyst, pH)

            switch catalyst
                case 'H2O'
                    concentration = 1000/18 * ones(size(pH));
                    pK_acceptor = obj.pK_catalysis('H3O+');
                case 'OH-'
                    concentration = 10.^(pH-14);
                    pK_acceptor = obj.pK_catalysis('H2O');
                case 'H2PO4'
                    [~, concentration, ~, ~] = obj.phosphateConcentration(pH);
                    pK_acceptor = obj.pK_catalysis('H3PO4');
                case 'HPO4'
                    [~, ~, concentration, ~] = obj.phosphateConcentration(pH);
                    pK_acceptor = obj.pK_catalysis('H2PO4');
                case 'PO4'
                    [~, ~, ~, concentration] = obj.phosphateConcentration(pH);
                    pK_acceptor = obj.pK_catalysis('HPO4');
            end

            krate = concentration .* obj.rateConstant(pK_donor, pK_acceptor);
            
        end
        
        % Reaction rate for reactions involving HA and H2A+
        function [krate] = rateReactionPerCatalyst2(obj, pK_acceptor, catalyst, pH)

            switch catalyst
                case 'H3O+'
                    concentration = 10.^(-pH);
                    pK_donor = obj.pK_catalysis('H3O+');
                case 'H2O'
                    concentration = 1000/18 * ones(size(pH));
                    pK_donor = obj.pK_catalysis('H2O');
                case 'H3PO4'
                    [concentration, ~, ~, ~] = obj.phosphateConcentration(pH);
                    pK_donor = obj.pK_catalysis('H3PO4');
                case 'H2PO4'
                    [~, concentration, ~, ~] = obj.phosphateConcentration(pH);
                    pK_donor = obj.pK_catalysis('H2PO4');
                case 'HPO4'
                    [~, ~, concentration, ~] = obj.phosphateConcentration(pH);
                    pK_donor = obj.pK_catalysis('HPO4');
            end

            krate = concentration .* obj.rateConstant(pK_donor, pK_acceptor);
            
        end
        
        % Reaction rate in a water solution for reactions involving HA and
        % A-
        function [krate] = rateReactionInWater(obj, pK_donor, pH)
            krate_H2O = obj.rateReactionPerCatalyst(pK_donor, 'H2O', pH);
            krate_OH = obj.rateReactionPerCatalyst(pK_donor, 'OH-', pH);
            krate = krate_H2O + krate_OH;
        end
        
        % Reaction rate in a PBS solution for reactions involving HA and
        % A-
        function [krate] = rateReactionInPBS(obj, pK_donor, pH)
            krate_water = obj.rateReactionInWater(pK_donor, pH);
            krate_H2PO4 = obj.rateReactionPerCatalyst(pK_donor, 'H2PO4', pH);
            krate_HPO4 = obj.rateReactionPerCatalyst(pK_donor, 'HPO4', pH);
            krate_PO4 = obj.rateReactionPerCatalyst(pK_donor, 'PO4', pH);

            krate = krate_water + krate_H2PO4 + krate_HPO4 + krate_PO4;
        end
        
        % Reaction rate in a water solution for reactions involving HA, A-,
        % and H2A+
        function [krate] = rateReactionInWater2(obj, pK_donor, pK_acceptor, pH)
            krate_deprotonation = obj.rateReactionInWater(pK_donor, pH);
            
            krate_H3O = obj.rateReactionPerCatalyst2(pK_acceptor, 'H3O+', pH);
            krate_H2O = obj.rateReactionPerCatalyst2(pK_acceptor, 'H2O', pH);
            
            krate = krate_deprotonation + krate_H3O + krate_H2O;
        end
        
        % Reaction rate in a PBS solution for reactions involving HA, A-,
        % and H2A+
        function [krate] = rateReactionInPBS2(obj, pK_donor, pK_acceptor, pH)
            krate_deprotonation = obj.rateReactionInPBS(pK_donor, pH);
            
            krate_H3O = obj.rateReactionPerCatalyst2(pK_acceptor, 'H3O+', pH);
            krate_H2O = obj.rateReactionPerCatalyst2(pK_acceptor, 'H2O', pH);
            krate_H3PO4 = obj.rateReactionPerCatalyst2(pK_acceptor, 'H3PO4', pH);
            krate_H2PO4 = obj.rateReactionPerCatalyst2(pK_acceptor, 'H2PO4', pH);
            krate_HPO4 = obj.rateReactionPerCatalyst2(pK_acceptor, 'HPO4', pH);

            krate = krate_deprotonation + krate_H3O + krate_H2O + krate_H3PO4 + krate_H2PO4 + krate_HPO4;
        end
        
        % Before calculating Z-spectra and MTRasym, set up the chemical
        % shift and relaxation times
        
        function obj = ParametersForBlochMcConnell(obj,chemical_shift, T1_I, T2_I, T1_S, T2_S)
            obj.chemical_shift = chemical_shift;
            obj.T1_I = T1_I;
            obj.T2_I = T2_I;
            obj.T1_S = T1_S;
            obj.T2_S = T2_S;
        end
        
        % Solving Bloch-McConnell equations
        %
        % Reference
        % Magnus Helgstrand, Torleif Härd & Peter Allard
        % J Biomol NMR. 2000;18:49?63.
        % DOI: 10.1023/A:1008309220156
        %
        function [PI, PS] = NumericalSolution(obj, freq_offsets, w1, sat_time, krate, concentration)
            
            [krate_, freq_] = meshgrid(krate, freq_offsets);
            freq_ = freq_(:);
            krate_ = krate_(:);
            
            concentration_ = repmat(concentration, [length(freq_offsets), 1]);
            concentration_ = concentration_(:);

            PI = zeros(1, length(krate_));
            PS = zeros(1, length(krate_));
            
            % Change 'parfor' to 'for' if you don't want to use parallel processing.
            parfor m=1:length(krate_)
                Omega = 2.0*pi*[freq_(m), freq_(m) - obj.chemical_shift];
                rho0 = [0.5 0 0 1 0 0 1]'; % Initial state
                
                P = [0, 0, 0, 0, 0, 0, 0;
                    0, 1/obj.T2_I + krate_(m) * concentration_(m), Omega(1), -w1(2), -krate_(m) * concentration_(m), 0, 0;
                    0, -Omega(1), 1/obj.T2_I + krate_(m) * concentration_(m), w1(1), 0, -krate_(m) * concentration_(m), 0;
                    -2.0/obj.T1_I, w1(2), -w1(1), 1/obj.T1_I + krate_(m) * concentration_(m), 0, 0, -krate_(m) * concentration_(m);
                    0, -krate_(m), 0, 0, 1/obj.T2_S + krate_(m), Omega(2), -w1(2);
                    0, 0, -krate_(m), 0, -Omega(2), 1/obj.T2_S + krate_(m), w1(1);
                    -2.0/obj.T1_S, 0, 0, -krate_(m), w1(2), -w1(1), 1/obj.T1_S + krate_(m)];
                
                rho1 = expm(-P * sat_time) * rho0;
                
                PI(1,m) = rho1(4);
                PS(1,m) = rho1(7);
            end
            
            PI = reshape(PI, [length(freq_offsets), length(krate)]);
            PS = reshape(PI, [length(freq_offsets), length(krate)]);
            
        end
        
        [concentration_] = HendersonHasselbach(obj, pK_donor, pH, concentration)

        [PI, PS] = ZSpec(obj, freq_offsets, w1, sat_time, pH, catalyst, pK_donor, concentration)
        
        [MTRasym, freq_offsets] = MTRasym(obj, freq_max, freq_step, w1, sat_time, pH, catalyst, pK_donor, concentration)
        
        [PI, PS] = ZSpec2(obj, freq_offsets, w1, sat_time, pH, catalyst, pK_donor, concentration, pK_acceptor)
        
        [MTRasym, freq_offsets] = MTRasym2(obj, freq_max, freq_step, w1, sat_time, pH, catalyst, pK_donor, concentration, pK_acceptor)
        
        [PI, PS] = ZSpecWater(obj, freq_offsets, w1, sat_time, pH, pK_donor, concentration)
        
        [PI, PS] = ZSpecWater2(obj, freq_offsets, w1, sat_time, pH, pK_donor, concentration, pK_acceptor)
        
        [MTRasym, freq_offsets] = MTRasymWater(obj, freq_max, freq_step, w1, sat_time, pH, pK_donor, concentration)
        
        [MTRasym, freq_offsets] = MTRasymWater2(obj, freq_max, freq_step, w1, sat_time, pH, pK_donor, concentration, pK_acceptor)
        
        [PI, PS] = ZSpecPBS(obj, freq_offsets, w1, sat_time, pH, pK_donor, concentration)
        
        [PI, PS] = ZSpecPBS2(obj, freq_offsets, w1, sat_time, pH, pK_donor, concentration, pK_acceptor)
        
        [MTRasym, freq_offsets] = MTRasymPBS(obj, freq_max, freq_step, w1, sat_time, pH, pK_donor, concentration)
        
        [MTRasym, freq_offsets] = MTRasymPBS2(obj, freq_max, freq_step, w1, sat_time, pH, pK_donor, concentration, pK_acceptor)
    end
    
end

