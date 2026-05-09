% Function to convert literature units to metal units
% Carter Reed

function metal = convertMetalUnits(metal)

Ry2eV  = 13.605693; % ev/Ry
nm2Ang = 10; % A/nm

% Convert SK
if isfield(metal, 'SK')
    metal.SK = convertStruct(metal.SK, Ry2eV);
end

% Convert onsite energies
if isfield(metal, 'onsite')
    metal.onsite = convertStruct(metal.onsite, Ry2eV);
end

% Convert lattice parameters to Angstroms
if isfield(metal, 'latParam')
    metal.latParam = convertStruct(metal.latParam, nm2Ang);
end

end


% Runs the actual unit conversion
function S = convertStruct(S, factor)
    fields = fieldnames(S);
    for i = 1:length(fields)
        f = fields{i};
        S.(f) = S.(f) * factor;
    end
end