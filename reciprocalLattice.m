% Function to compute reciprocal lattice vectors
% Carter Reed

function b = reciprocalLattice(a)

% For stability, force double precision
a = double(a);

% a is a 3x3 matrix where each row is a primitive lattice vector
a1 = a(1,:);
a2 = a(2,:);
a3 = a(3,:);

% Unit cell volume
V = dot(a1, cross(a2, a3));

% Reciprocal lattice vectors
b1 = 2*pi*cross(a2, a3)/V;
b2 = 2*pi*cross(a3, a1)/V;
b3 = 2*pi*cross(a1, a2)/V;

% Compress into 3x3 matrix
b = [b1; b2; b3];

end