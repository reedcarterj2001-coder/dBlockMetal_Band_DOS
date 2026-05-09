% Function to build eigenvalue cache for optimization
% Carter Reed

function spectrum = computeSpectrum(crystal, kvec)

Nk = size(kvec, 1);
H0 = buildHamiltonian(crystal, kvec(1,:));
nBands = size(H0,1);

bands = zeros(nBands, Nk);
eigvecs = zeros(nBands, nBands, Nk);

allE = zeros(Nk*nBands,1);
idx = 1;

for i = 1:Nk

    H = buildHamiltonian(crystal, kvec(i,:));
    [V, D] = eig(H);
    E = diag(D);

    % Sort consistently
    [E, order] = sort(real(E));
    V = V(:,order);

    bands(:,i) = E;
    eigvecs(:,:,i) = V;

    allE(idx:idx+nBands-1) = E;
    idx = idx + nBands;

end

% Package nicely
spectrum.bands   = bands;
spectrum.eigvecs = eigvecs;
spectrum.allE    = allE;
spectrum.Nk      = Nk;

end