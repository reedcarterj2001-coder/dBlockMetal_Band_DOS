% Function to plot band structure and DOS
% Carter Reed

function plotBandDOS(kdist, spectrum, Egrid, DOS, EF, tickPos, ...
                     tickLabels, crystal, sigma)

% Shift energies to make EF = 0 reference
bands = spectrum.bands - EF;
Eplot = Egrid - EF;

% Orbital wieghts
W = computeOrbitalWeights(spectrum);

% Number of bands
nBands = size(bands,1);
Nk = size(bands,2);

% -----------
% Make figure
% -----------

figure

str = [crystal.name, ' (', crystal.structure,')'];
sgtitle(str)

ylimVals = [min(Eplot), max(Eplot)];


% --------------
% Band Structure
% --------------

subplot(1,2,1)
hold on
for n = 1:nBands
    for i = 1:Nk-1

        % Normalized weights
        ws = W.s(n,i);
        wp = W.p(n,i);
        wd = W.d(n,i);

        % Color of band where
        % pure s = red
        % pure p = green
        % pure d = blue
        color = ws*[1 0 0] + wp*[0 1 0] + wd*[0 0 1];

        plot(kdist(i:i+1), bands(n,i:i+1), 'LineWidth', 2, 'Color', color);

    end
end

% Plot Fermi level
yline(0, 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');

% Vertical symmetry lines
for i = 1:length(tickPos)
    xline(tickPos(i), '--', 'Color', [0.5 0.5 0.5]);
end

% Axis formatting
xlim([kdist(1), kdist(end)])
xticks(tickPos)
tickLabels = strrep(tickLabels,'Gamma','\Gamma');
xticklabels(tickLabels)

ylabel('Energy - E_F (eV)')

title('Fat-Band Structure')

ylim(ylimVals);

hold off


% -----------------
% Density of States
% -----------------

subplot(1,2,2)
hold on


DOS   = DOS(:);
Eplot = Eplot(:);

% Plot DOS
plot(DOS, Eplot, 'k', 'LineWidth', 2)

% Plot Fermi level
yline(0, 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-.');

% Shade occupied states
mask = (Eplot <= 0);
if any(mask)

    xfill = [zeros(sum(mask),1); flipud(DOS(mask))];
    yfill = [Eplot(mask); flipud(Eplot(mask))];

    fill(xfill, yfill, [0.2 0.6 1.0], ...
        'FaceAlpha',0.4, ...
        'EdgeColor','none');

end

xlabel('DOS')
ylabel('Energy - E_F (eV)')

str1 = ['Density of States with Gaussian Smearing (\sigma = ', ...
        num2str(sigma),' eV)']; 
title(str1)

ylim(ylimVals);
xlim auto

outdir = ; % Your desried figure directory here, put ' ' to leave blank
filename = fullfile(outdir,[crystal.name,'_', ...
                            crystal.structure,'_BandPlot_DOSPlot']);
exportgraphics(gcf, filename + ".tif",'Resolution', 300);