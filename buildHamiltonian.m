% Function to build system 6x6 Slater-Koster Hamiltonian in tight binding 
% model consisting of the s and d orbitals for a total of 6 orbitals
% Carter Reed

function H = buildHamiltonian(crystal,k)

% Orbitals: s, px, py, pz, dxy, dyz, dxz, dx^2-y^2, dz^2
norb = 9;

% Initialize Hamiltonian
H = zeros(norb, norb);

% -------------
% On-site terms
% -------------
O = crystal.onsite;
structure = crystal.structure;
H(1,1) = O.s;

% Select appropriate splitting for lattice type
switch upper(structure)
    case 'HCP'
        H(2,2) = O.p;
        H(3,3) = O.p;
        H(4,4) = O.p;
        H(5,5) = O.d1;
        H(6,6) = O.d1;
        H(7,7) = O.d1;
        H(8,8) = O.d2;
        H(9,9) = O.d0;

    case {'FCC','BCC'}
        H(2,2) = O.p;
        H(3,3) = O.p;
        H(4,4) = O.p;
        H(5,5) = O.d1;
        H(6,6) = O.d1;
        H(7,7) = O.d1;
        H(8,8) = O.d2;
        H(9,9) = O.d2;

    otherwise
        error('Unknown crystal structure: %s',structure);
end

% Reciprocal lattice vectors
b = crystal.reciprocalVectors;

% Cartesian wavevector
kcart = k(1)*b(1,:) + k(2)*b(2,:) + k(3)*b(3,:);

% -------------------------------
% Loop over neighbors for hopping
% -------------------------------
params = crystal.SK;
for n = 1:size(crystal.neighbors,1)
    rvec = crystal.neighbors(n,:);
    r = norm(rvec);
    rhat = rvec/r;
    lx = rhat(1);
    my = rhat(2);
    nz = rhat(3);
    phase = exp(1i * dot(kcart,rvec));

    % Radial factor
    r0 = crystal.r0;
    rc = 1.0;
    scale = exp(-(r - r0)/rc);

    T = SlaterKosterMatrix(lx, my, nz, params, norb);
    T = scale*T;
    H = H + T*phase;

end

% Hermitian clean up
H = (H + H')/2;

end
