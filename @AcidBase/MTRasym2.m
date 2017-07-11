function [MTRasym, freq_offsets] = MTRasym2(obj, freq_max, freq_step, w1, sat_time, pH, catalyst, pK_donor, concentration, pK_acceptor)

freq_offsets = -freq_max:freq_step:freq_max;
PI = obj.ZSpec2(freq_offsets, w1, sat_time, pH, catalyst, pK_donor, concentration, pK_acceptor);

MTRasym = flipud(PI) - PI;
end