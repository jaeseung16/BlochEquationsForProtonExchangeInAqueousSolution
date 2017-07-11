function [PI, PS] = ZSpec(obj, freq_offsets, w1, sat_time, pH, catalyst, pK_donor, concentration)
% w1 = [w1x, w1y] in rad/s

concentration_ = obj.HendersonHasselbach(pK_donor, pH, concentration);

% exchange rate
krate = obj.rateReactionPerCatalyst(pK_donor, catalyst, pH);

[PI, PS] = obj.NumericalSolution(freq_offsets, w1, sat_time, krate, concentration_);

end
