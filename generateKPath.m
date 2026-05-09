% Function to generate k-path for the coming band structure calculations
% Carter Reed

function [kvec, kdist, tickPos, tickLabels] = generateKPath(crystal, nPerSegment)

% Number of k-points between high symmetry points
if (nargin < 2)
    nPerSegment = 50;
end

path = crystal.kpath;
kpoints = crystal.kpoints;

kvec = [];
kdist = [];
tickPos = [];
tickLabels = {};

% Distance traveled along kpath
currentDist = 0;

for i = 1:length(path)-1

    k1 = kpoints.(path{i});
    k2 = kpoints.(path{i+1});

    % Interpolate linearly in k-space
    seg = zeros(nPerSegment,3);
    for j = 1:nPerSegment
        t = (j-1)/(nPerSegment-1);
        seg(j,:) = (1-t)*k1 + t*k2;
    end
    
    % Remove duplicate point between segemnts
    if (i > 1)
        seg = seg(2:end,:);
    end

    % Build distance for segment
    segDist = zeros(size(seg, 1), 1);

    % Compute distance
    for j = 2:size(seg,1)
        dk = seg(j,:) - seg(j-1,:);
        currentDist = currentDist + norm(dk);
        segDist(j) = currentDist;
    end

    % First point correction
    if isempty(kdist)
        segDist(1) = 0;
    else
        segDist(1) = kdist(end);
    end

    % Append consistently
    kvec  = [kvec ; seg];
    kdist = [kdist; segDist];

    % Store tick positions (start of each segment)
    tickPos = [tickPos; kdist(end)];
    tickLabels = [tickLabels; path{i}];

end

% Add final label
tickLabels = [tickLabels; path{end}];

end