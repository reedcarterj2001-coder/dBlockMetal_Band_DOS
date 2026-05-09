% Function to calculate orbital weights
% Carter Reed

function W = computeOrbitalWeights(spectrum)

eigvecs = spectrum.eigvecs;

% Indices
idx.s = 1;
idx.p = 2:4;
idx.d = 5:9;

[norb, nBands,Nk] = size(eigvecs);

% Initialize weights
W.s = zeros(nBands, Nk);
W.p = zeros(nBands, Nk);
W.d = zeros(nBands, Nk);

for ik = 1:Nk
    for ib = 1:nBands

        psi = spectrum.eigvecs(:,ib,ik); 

        % Raw weight
        ws = sum(abs(psi(idx.s)).^2);
        wp = sum(abs(psi(idx.p)).^2);
        wd = sum(abs(psi(idx.d)).^2);

        % Normalization factor
        wnorm = ws + wp + wd;

        % Normalized weights
        W.s(ib,ik) = ws / wnorm;
        W.p(ib,ik) = wp / wnorm;
        W.d(ib,ik) = wd / wnorm;

    end
end

end