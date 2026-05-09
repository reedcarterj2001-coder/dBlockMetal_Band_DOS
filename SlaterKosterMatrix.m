% Function to drive the Slater-Koster matrix elements
% Carter Reed

function T = SlaterKosterMatrix(l, m, n, p, norb)

% ------------------------------------------------------------------------
% Indices, Parameters, and Preallocation
% --------------------------------------
idx.s     = 1;
idx.px    = 2;
idx.py    = 3;
idx.pz    = 4;
idx.dxy   = 5;
idx.dyz   = 6;
idx.dxz   = 7;
idx.dx2y2 = 8;
idx.dz2   = 9;

% p = parameter structure
Vss    = p.Vsssig;
Vsp    = p.Vspsig;
Vppsig = p.Vppsig;
Vpppi  = p.Vpppi;
Vpdsig = p.Vpdsig;
Vpdpi  = p.Vpdpi;
Vsd    = p.Vsdsig;
Vddsig = p.Vddsig;
Vddpi  = p.Vddpi;
Vdddel = p.Vdddel;

T = zeros(norb,norb);

% ------------------------------------------------------------------------
% s-s coupling
% ------------
T(idx.s,idx.s) = Vss;

% ------------------------------------------------------------------------
% s-p coupling
% ------------
T(idx.s,idx.px) = l*Vsp; % px
T(idx.s,idx.py) = m*Vsp; % py
T(idx.s,idx.pz) = n*Vsp; % pz

% Hermitian symmetry
T(idx.px:idx.pz,idx.s) = T(idx.s,idx.px:idx.pz);

% ------------------------------------------------------------------------
% p-p coupling
% ------------
T(idx.px,idx.px) = Vppsig*(l^2) + Vpppi*(1-l^2);
T(idx.py,idx.py) = Vppsig*(m^2) + Vpppi*(1-m^2);
T(idx.pz,idx.pz) = Vppsig*(n^2) + Vpppi*(1-n^2);

T(idx.px,idx.py) = (Vppsig - Vpppi)*l*m;
T(idx.px,idx.pz) = (Vppsig - Vpppi)*l*n;
T(idx.py,idx.pz) = (Vppsig - Vpppi)*m*n;

% Hermitian symmetry
T(idx.py,idx.px) = T(idx.px,idx.py);
T(idx.pz,idx.px) = T(idx.px,idx.pz);
T(idx.pz,idx.py) = T(idx.py,idx.pz);

% ------------------------------------------------------------------------
% s-d coupling
% ------------
T(idx.s,idx.dxy)   = Vsd*sqrt(3)*l*m;   %dxy
T(idx.s,idx.dyz)   = Vsd*sqrt(3)*m*n;   %dyz
T(idx.s,idx.dxz)   = Vsd*sqrt(3)*l*n;   %dxz
T(idx.s,idx.dx2y2) = Vsd*(l^2 - m^2)/2; %dx^2-y^2
T(idx.s,idx.dz2)   = Vsd*(3*n^2 - 1)/2; %dz^2

% Hermitian symmetry
T(idx.dxy:idx.dz2,idx.s) = T(idx.s,idx.dxy:idx.dz2);

% ------------------------------------------------------------------------
% p-d coupling
% ------------

% px coupled to
% -------------

% dxy
T(idx.px,idx.dxy)   = ...
            Vpdsig*sqrt(3)*(l^2)*m           + Vpdpi*m*(1 - 2*l^2);
% dyz
T(idx.px,idx.dyz)   = ...
            Vpdsig*sqrt(3)*l*m*n             - Vpdpi*2*l*m*n;
% dxz
T(idx.px,idx.dxz)   = ...
            Vpdsig*sqrt(3)*(l^2)*n           + Vpdpi*n*(1 - 2*l^2);
% dx2y2
T(idx.px,idx.dx2y2) = ...
            Vpdsig*(sqrt(3)/2)*l*(l^2 - m^2) + Vpdpi*l*(1 - l^2 + m^2);
% dz2
T(idx.px,idx.dz2)   = ...
            Vpdsig*l*(n^2 - 0.5*(l^2 + m^2)) - Vpdpi*sqrt(3)*l*(n^2);

% Hermitian symmetry
T(idx.dxy:idx.dz2,idx.px) = T(idx.px,idx.dxy:idx.dz2);


% py coupled to
% -------------

% dxy
T(idx.py,idx.dxy)   = ...
            Vpdsig*sqrt(3)*(m^2)*l           + Vpdpi*l*(1 - 2*m^2);
% dyz
T(idx.py,idx.dyz)   = ...
            Vpdsig*sqrt(3)*(m^2)*n           + Vpdpi*n*(1 - 2*m^2);
% dxz
T(idx.py,idx.dxz)   = ...
            Vpdsig*sqrt(3)*l*m*n             - Vpdpi*2*l*m*n;
% dx2y2
T(idx.py,idx.dx2y2) = ...
            Vpdsig*(sqrt(3)/2)*m*(m^2 - l^2) + Vpdpi*m*(1 + l^2 - m^2);
% dz2
T(idx.py,idx.dz2)   = ...
            Vpdsig*m*(n^2 - 0.5*(l^2 + m^2)) - Vpdpi*sqrt(3)*m*(n^2);

% Hermitian symmetry
T(idx.dxy:idx.dz2,idx.py) = T(idx.py,idx.dxy:idx.dz2);


% pz coupled to
% -------------

% dxy
T(idx.pz,idx.dxy)   = ...
            Vpdsig*sqrt(3)*l*m*n             + Vpdpi*n*(1 - 2*l^2);
% dyz
T(idx.pz,idx.dyz)   = ...
            Vpdsig*sqrt(3)*m*(n^2)           + Vpdpi*m*(1 - 2*n^2);
% dxz
T(idx.pz,idx.dxz)   = ...
            Vpdsig*sqrt(3)*l*(n^2)           + Vpdpi*l*(1 - 2*n^2);
% dx2y2
T(idx.pz,idx.dx2y2) = ...
            Vpdsig*(sqrt(3)/2)*n*(l^2 - m^2) - Vpdpi*n*(l^2 - m^2);
% dz2
T(idx.pz,idx.dz2)   = ...
            Vpdsig*n*(n^2 - 0.5*(l^2 + m^2)) + Vpdpi*sqrt(3)*n*(l^2 + m^2);

% Hermitian symmetry
T(idx.dxy:idx.dz2,idx.pz) = T(idx.pz,idx.dxy:idx.dz2);


% ------------------------------------------------------------------------
% d-d coupling
% ------------
% dxy - dxy
T(idx.dxy,idx.dxy) = ...
            Vddsig*(l^2)*(m^2)   + Vddpi*(l^2 + m^2) + Vdddel*(n^2);
% dyz - dyz
T(idx.dyz,idx.dyz) = ...
            Vddsig*(m^2)*(n^2)   + Vddpi*(m^2 + n^2) + Vdddel*(l^2);
% dxz - dxz
T(idx.dxz,idx.dxz) = ...
            Vddsig*(l^2)*(n^2)   + Vddpi*(l^2 + n^2) + Vdddel*(m^2);
% dx2y2 - dx2y2
T(idx.dx2y2,idx.dx2y2) = ...
            Vddsig*(l^2 - m^2)^2 + Vddpi*(l^2 + m^2) + Vdddel*(n^2);
% dz2 - dz2
T(idx.dz2,idx.dz2) = ...
            Vddsig*(3*n^2 - 1)   + Vddpi*(1 - n^2)   + Vdddel*(n^2);

T(idx.dxy,idx.dyz) = Vddsig* l * (m^2)*n;    % dxy - dyz
T(idx.dxy,idx.dxz) = Vddsig*(l^2)*m *  n;    % dxy - dxz
T(idx.dyz,idx.dxz) = Vddsig* l *  m * (n^2); % dyz - dxz

% Hermitian symmetry
T(idx.dyz,idx.dxy) = T(idx.dxy,idx.dyz);
T(idx.dxz,idx.dxy) = T(idx.dxy,idx.dxz);
T(idx.dxz,idx.dyz) = T(idx.dyz,idx.dxz);

% ------------------------------------------------------------------------
% Keep remaining off-diagonal d-d terms minimal
T = (T + T')/2;

end