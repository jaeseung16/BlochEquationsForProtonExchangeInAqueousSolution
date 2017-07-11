function [MTRasym, freq_offsets] = MTRasymWater(obj, freq_max, freq_step, w1, sat_time, pH, pK_donor, concentration)
freq_offsets = -freq_max:freq_step:freq_max;

[PI, PS] = obj.ZSpecWater(freq_offsets, w1, sat_time, pH, pK_donor, concentration);

MTRasym = flipud(PI) - PI;
end