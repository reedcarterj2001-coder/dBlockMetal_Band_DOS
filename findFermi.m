% Function to find Fermi level
% Carter Reed

function EF = findFermi(Egrid, DOS, crystal)

% Number of electrons per primitive cell
Ne = crystal.Valence * crystal.nAtoms;

Egrid = Egrid(:);
DOS   = DOS(:);

% Ensure sorted grid
[Egrid, idx] = sort(Egrid);
DOS = DOS(idx);

% Cumulative electron count
Ncum = cumtrapz(Egrid, DOS);

% Force strict monotonicity
Ncum = Ncum + 1e-12*(1:length(Ncum))';

% Remove duplicates safely
[Ncum_unique, ia] = unique(Ncum, 'stable');
E_unique = Egrid(ia);

% Invert
EF = interp1(Ncum_unique, E_unique, Ne, 'linear', 'extrap');

end