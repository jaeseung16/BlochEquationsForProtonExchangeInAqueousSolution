function [rateF] = rateConstant(obj, pK_donor, pK_acceptor)
rateF = obj.k0 / ( 1 + 10^(pK_donor - pK_acceptor) );
end
