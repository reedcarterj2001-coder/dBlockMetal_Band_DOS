% Function to build the crystal systems
% High symmetry k-points taken from Wikipedia
% (https://en.wikipedia.org/wiki/Brillouin_zone)
% Other crystal structures may be added in future updates.
% Crystal systems using primitive reciprocal basis.
% Carter Reed

function crystal = buildCrystal(metal)

% Copy metal data
crystal = metal;

% Select structure type
structure = upper(metal.structure);

% Select appropriate structure
switch structure

    case 'BCC'
        a = metal.latParam.a;

        % Primitive lattice vectors
        crystal.realVectors = (a/2)*[-1,  1,  1;
                                      1, -1,  1;
                                      1,  1, -1];

        % Crystal Basis (one atome at origin)
        crystal.primitiveBasis = [0, 0, 0];

        % High symmetry k-points
        crystal.kpoints.Gamma = [0   , 0   , 0   ];
        crystal.kpoints.H     = [0   , 0   , 1   ];
        crystal.kpoints.N     = [0   , 0.5 , 0   ];
        crystal.kpoints.P     = [0.25, 0.25, 0.25];

        crystal.kpath = {'Gamma','H','N','Gamma','P','H'};

        % Nearest Neighbors
        crystal.neighbors = (a/2)*[ 1,  1,  1;
                                   -1,  1,  1;
                                    1, -1,  1;
                                    1,  1, -1;
                                    1, -1, -1;
                                   -1,  1, -1;
                                   -1, -1,  1;
                                   -1, -1, -1];

        % Nearest neighbor distance
        crystal.r0 = min(vecnorm(crystal.neighbors,2,2));

    case 'FCC'
        a = metal.latParam.a;

        % Primitive lattice vectors
        crystal.realVectors = (a/2)*[0, 1, 1;
                                     1, 0, 1;
                                     1, 1, 0];

        % Crystal Basis (one atom at origin)
        crystal.primitiveBasis = [0, 0, 0];

        % High symmetry k-points
        crystal.kpoints.Gamma = [0    , 0   , 0    ];
        crystal.kpoints.L     = [0.5  , 0.5 , 0.5  ];
        crystal.kpoints.K     = [0.375, 0.75, 0.375];
        crystal.kpoints.W     = [1    , 0.5 , 0    ];
        crystal.kpoints.U     = [1    , 0.25, 0.25 ];
        crystal.kpoints.X     = [1    , 0   , 0    ];

        crystal.kpath = {'Gamma','X','W','K','Gamma','L','U','X'};

        % Nearest Neighbors
        crystal.neighbors = (a/2)*[ 0,  1,  1;
                                    0, -1,  1;
                                    0,  1, -1;
                                    0, -1, -1;
                                    1,  0,  1;
                                   -1,  0,  1;
                                    1,  0, -1;
                                   -1,  0, -1;
                                    1,  1,  0;
                                   -1,  1,  0;
                                    1, -1,  0;
                                   -1, -1,  0];

        % Nearest neighbor distance
        crystal.r0 = min(vecnorm(crystal.neighbors,2,2));

    case 'HCP'
        a = metal.latParam.a;
        c = metal.latParam.c;

        % Primitive lattice vectors
        crystal.realVectors = [a   , 0          , 0;
                               -a/2, sqrt(3)*a/2, 0;
                               0   , 0          , c];

        % Crystal Basis (two-atom basis)
        crystal.primitiveBasis = [0  , 0  , 0  ;
                                  2/3, 1/3, 1/2];

        % High Symmetry k-points
        crystal.kpoints.Gamma = [0  , 0  , 0  ];
        crystal.kpoints.M     = [0.5, 0  , 0  ];
        crystal.kpoints.K     = [1/3, 1/3, 0  ];
        crystal.kpoints.A     = [0  , 0  , 0.5];

        crystal.kpath = {'Gamma','M','K','Gamma','A'};

        % Nearest Neighbors
        r1  = crystal.realVectors(1,:);
        r2  = -r1;
        r3  = crystal.realVectors(2,:);
        r4  = -r3;
        r5  = crystal.realVectors(1,:)+crystal.realVectors(2,:);
        r6  = -r5;
        r7  =  (2/3)*r1 + (1/3)*r3 + (1/2)*r5;
        r8  = -r7;
        r9  = -(1/3)*r1 + (1/3)*r3 + (1/2)*r5;
        r10 = -r9;
        r11 =  (2/3)*r1 - (2/3)*r3 + (1/2)*r5;
        r12 = -r11;
        crystal.neighbors = [r1 ;
                             r2 ;
                             r3 ;
                             r4 ;
                             r5 ;
                             r6 ;
                             r7 ;
                             r8 ;
                             r9 ;
                             r10;
                             r11;
                             r12];

        % Nearest neighbor distance
        crystal.r0 = min(vecnorm(crystal.neighbors,2,2));

    otherwise
        error("Unknown crystal structure: %s",metal.structure);

end

% Compute number of atoms per primitive cell
crystal.nAtoms = size(crystal.primitiveBasis,1);

% Compute reciprocal lattice vectors
crystal.reciprocalVectors = reciprocalLattice(crystal.realVectors);

end