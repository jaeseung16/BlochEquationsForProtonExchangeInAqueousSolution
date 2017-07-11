function [PI, PS] = ZSpecPBS2(obj, freq_offsets, w1, sat_time, pH, pK_donor, concentration, pK_acceptor)

% w1 = [w1x, w1y] in rad/s

concentration_ = obj.HendersonHasselbach(pK_donor, pH, concentration);

% exchange rate
krate = obj.rateReactionInPBS2(pK_donor, pK_acceptor, pH);

[PI, PS] = obj.NumericalSolution(freq_offsets, w1, sat_time, krate, concentration_);

end
