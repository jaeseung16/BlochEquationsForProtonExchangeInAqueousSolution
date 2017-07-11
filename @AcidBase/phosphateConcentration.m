function [phos1, phos2, phos3, phos4] = phosphateConcentration(obj, pH)
K = 10.^pH;

factor = 1 + (obj.K1) .* K .* (1 + (obj.K2) .* K .* (1 + (obj.K3) .* K ) );

phos1 = (obj.phosphate_concentration) ./ factor; % [H3PO4]
phos2 = phos1 .* obj.K1 .* K; % [H2PO4-]
phos3 = phos2 .* obj.K2 .* K; % [HPO42-]
phos4 = phos3 .* obj.K3 .* K; % [HPO43-]
end