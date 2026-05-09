% Function to compute DOS
% Carter Reed

function [Egrid, DOS] = computeDOS(spectrum, nBins, sigma)

allE = spectrum.allE(:);
Nk = spectrum.Nk;

pad = 5*sigma;
Emin = min(allE) - pad;
Emax = max(allE) + pad;

Egrid = linspace(Emin, Emax, nBins);
DOS = zeros(size(Egrid));

A = (1/(sqrt(2*pi)*sigma)); % Gaussian normalization constant
for i = 1:length(allE)
    DOS = DOS + A*exp(-(Egrid - allE(i)).^2 / (2*sigma^2));
end

% Account for spin and normalize
DOS = 2*DOS/Nk;

DOS = DOS(:);

end