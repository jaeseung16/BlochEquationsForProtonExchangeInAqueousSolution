function [MTRasym, freq_offsets] = MTRasym(obj, freq_max, freq_step, w1, sat_time, pH, catalyst, pK_donor, concentration)

freq_offsets = -freq_max:freq_step:freq_max;

[PI, PS] = obj.ZSpec(freq_offsets, w1, sat_time, pH, catalyst, pK_donor, concentration);

MTRasym = flipud(PI) - PI;

end
