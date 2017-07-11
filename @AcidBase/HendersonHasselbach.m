function [concentration_] = HendersonHasselbach(obj, pK_donor, pH, concentration)
concentration_HH = concentration ./ (1 + 10.^(pH - pK_donor) );
concentration_ = concentration_HH / (1000/18); % [HA] / [Water]
end