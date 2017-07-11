function [PI, PS] = ZSpec2(obj, freq_offsets, w1, sat_time, pH, catalyst, pK_donor, concentration, pK_acceptor)
% w1 = [w1x, w1y] in rad/s

concentration_ = obj.HendersonHasselbach(pK_donor, pH, concentration);

% exchange rate
krate = obj.rateReactionPerCatalyst2(pK_acceptor, catalyst, pH);

[PI, PS] = obj.NumericalSolution(freq_offsets, w1, sat_time, krate, concentration_);

end
