% Function to pull information from metalDatabase.m
% Carter Reed

function metal = getTransitionMetal(symbol)

% Load metalDatabase.m
db = metalDatabase();

% Force user input to lowercase
symbol = lower(symbol);

% Vaildate entry
if isfield(db, symbol)
    metal = db.(symbol);
else
    error("d-Block Transition Metal not found.");
end

end