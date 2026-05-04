% NACA 4-digit wind-tunnel simulator using a doublet panel method.

clear
clc
close all

fprintf('NACA 4-digit wind-tunnel simulator\n');
NACAcode = strtrim(input('Enter the NACA 4 digit code: ', 's'));
validateNacaCode(NACAcode);

if strcmp(NACAcode, '2412')
    runNaca2412Study(NACAcode);
else
    uinf = readPositiveNumber('Enter the freestream velocity in m/s: ');
    alpha = readNumber('Enter the angle of attack in degrees: ');
    N = readPositiveInteger('Enter the number of panels: ');

    result = solvePanelMethod(NACAcode, alpha, N, uinf);
    fprintf('The coefficient of lift is %.4f\n', result.Cl);

    plotFlowField(result, 300, 'Streamline plot.png', 'Arrow plot.png');
end

function runNaca2412Study(NACAcode)
uinf = 15;
panelCounts = [50, 100, 200];
alphas = 0:2:10;
CL = zeros(numel(panelCounts), numel(alphas));

for i = 1:numel(panelCounts)
    for j = 1:numel(alphas)
        result = solvePanelMethod(NACAcode, alphas(j), panelCounts(i), uinf);
        CL(i, j) = result.Cl;
    end
end

plotLiftCurves(alphas, CL, panelCounts);

flowResult = solvePanelMethod(NACAcode, 10, 200, uinf);
plotFlowField(flowResult, 200, ...
    'Streamline plot NACA-2412.png', ...
    'Arrow plot for NACA-2412.png');
end

function result = solvePanelMethod(NACAcode, alpha, N, uinf)
[X, Z] = panelgen(NACAcode, alpha, N);
beta = atan2(diff(Z(1:N+1)), diff(X(1:N+1)));
A = zeros(N + 1, N + 1);
b = zeros(N + 1, 1);

for i = 1:N
    p = [(X(i + 1) + X(i)) / 2, (Z(i + 1) + Z(i)) / 2];
    b(i) = -uinf * sin(deg2rad(alpha) - beta(i));

    for j = 1:N+1
        p1 = [X(j), Z(j)];
        p2 = [X(j + 1), Z(j + 1)];
        [U, V] = cdoublet(p, p1, p2);
        A(i, j) = V * cos(beta(i)) - U * sin(beta(i));
    end
end

% Kutta condition.
A(N + 1, 1) = 1;
A(N + 1, N) = -1;
A(N + 1, N + 1) = 1;

mu = A \ b;

result.NACAcode = NACAcode;
result.alpha = alpha;
result.N = N;
result.uinf = uinf;
result.X = X;
result.Z = Z;
result.mu = mu;
result.Cl = -2 * mu(end) / uinf;
end

function plotLiftCurves(alphas, CL, panelCounts)
xfoilFile = 'xf-naca2412-il-1000000.txt';
xfoilData = readmatrix(xfoilFile, 'FileType', 'text');
validRows = all(isfinite(xfoilData(:, 1:2)), 2);
xfoilData = xfoilData(validRows, 1:2);

for i = 1:numel(panelCounts)
    figure(i)
    plot(xfoilData(:, 1), xfoilData(:, 2), 'k', 'DisplayName', 'XFOIL data')
    hold on
    plot(alphas, CL(i, :), 'o-', ...
        'DisplayName', sprintf('Panel method (%d panels)', panelCounts(i)))
    hold off
    grid on
    xlim([0, 10])
    title(sprintf('Graph of Cl vs AoA for %d panels', panelCounts(i)))
    xlabel('Angle of attack (degrees)')
    ylabel('Coefficient of lift')
    legend('Location', 'best')
    saveas(gcf, sprintf('Cl vs AoA - %d panels.png', panelCounts(i)))
end
end

function plotFlowField(result, gridSize, streamlineFile, arrowFile)
xGrid = linspace(-0.2, 1.2, gridSize);
zGrid = linspace(-0.7, 0.7, gridSize);
uGrid = zeros(gridSize);
vGrid = zeros(gridSize);

for row = 1:gridSize
    for col = 1:gridSize
        p = [xGrid(col), zGrid(row)];
        uGrid(row, col) = result.uinf * cosd(result.alpha);
        vGrid(row, col) = result.uinf * sind(result.alpha);

        for panel = 1:result.N+1
            p1 = [result.X(panel), result.Z(panel)];
            p2 = [result.X(panel + 1), result.Z(panel + 1)];
            [U, V] = cdoublet(p, p1, p2);
            uGrid(row, col) = uGrid(row, col) + result.mu(panel) * U;
            vGrid(row, col) = vGrid(row, col) + result.mu(panel) * V;
        end
    end
end

[Xg, Zg] = meshgrid(xGrid, zGrid);
insideAirfoil = inpolygon(Xg, Zg, result.X(1:result.N), result.Z(1:result.N));
Xg(insideAirfoil) = NaN;
Zg(insideAirfoil) = NaN;
uGrid(insideAirfoil) = NaN;
vGrid(insideAirfoil) = NaN;

figure
plot(result.X(1:end-1), result.Z(1:end-1), 'LineWidth', 1.5)
hold on
streamslice(Xg, Zg, uGrid, vGrid)
hold off
formatFlowPlot(result)
title(sprintf('Streamline plot for NACA-%s at %.1f degrees with %d panels', ...
    result.NACAcode, result.alpha, result.N))
saveas(gcf, streamlineFile)

sample = 1:20:gridSize;
figure
plot(result.X(1:end-1), result.Z(1:end-1), 'LineWidth', 1.5)
hold on
quiver(Xg(sample, sample), Zg(sample, sample), ...
    uGrid(sample, sample), vGrid(sample, sample))
hold off
formatFlowPlot(result)
title(sprintf('Arrow plot for NACA-%s at %.1f degrees with %d panels', ...
    result.NACAcode, result.alpha, result.N))
saveas(gcf, arrowFile)
end

function formatFlowPlot(result)
grid on
axis equal
xlim([-0.2, 1.2])
ylim([-0.7, 0.7])
xlabel('x/c')
ylabel('z/c')
text(-0.18, 0.62, sprintf('U = %.1f m/s, Cl = %.4f', result.uinf, result.Cl))
end

function validateNacaCode(NACAcode)
if numel(NACAcode) ~= 4 || ~all(isstrprop(NACAcode, 'digit'))
    error('Please enter a valid four digit NACA code, for example 2412.');
end
end

function value = readNumber(prompt)
value = input(prompt);
if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value)
    error('Input must be a finite number.');
end
end

function value = readPositiveNumber(prompt)
value = readNumber(prompt);
if value <= 0
    error('Input must be greater than zero.');
end
end

function value = readPositiveInteger(prompt)
value = readPositiveNumber(prompt);
if value ~= floor(value)
    error('Input must be a whole number.');
end
end
