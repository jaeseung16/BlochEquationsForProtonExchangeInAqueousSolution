function [MTRasym, freq_offsets] = MTRasymPBS2(obj, freq_max, freq_step, w1, sat_time, pH, pK_donor, concentration, pK_acceptor)

freq_offsets = -freq_max:freq_step:freq_max;

PI = obj.ZSpecPBS2(freq_offsets, w1, sat_time, pH, pK_donor, concentration, pK_acceptor);

MTRasym = flipud(PI) - PI;

end
