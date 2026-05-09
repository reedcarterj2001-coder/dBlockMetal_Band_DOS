% User input script for band structure calculations
% Carter Reed

clc
clear

% Prompt user for the atomic symbol of their desired d-block transition
% metal
symbol = input("Please enter the atomic symbol of your desired d-block transition metal: ","s");

% Retrieve data
metal = getTransitionMetal(symbol);
% Convert to metal units
metal = convertMetalUnits(metal);

% Build crystal system
crystal = buildCrystal(metal);

% Generate k-path for selected system
nPerSegment = 100; % Number of k-points between high symmetry points
[kvec, kdist, tickPos, tickLabels] = generateKPath(crystal, nPerSegment);

% Compute spectrum
spectrum = computeSpectrum(crystal, kvec);

% Compute DOS
nBins = 3000; % How finely sample the E axis is for DOS
sigma = 0.07;  % Gaussian smearing parameter in eV
[Egrid, DOS] = computeDOS(spectrum, nBins, sigma);

% Find Fermi level
EF = findFermi(Egrid, DOS, crystal);

% Plot the bands and DOS
plotBandDOS(kdist, spectrum, Egrid, DOS, EF, tickPos, tickLabels, ...
            crystal, sigma);