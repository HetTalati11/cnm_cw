function [X, Z] = panelgen(NACAcode, alpha, N)
%PANELGEN Discretise a NACA 4-digit aerofoil using cosine-spaced panels.
%   [X, Z] = PANELGEN(NACAcode, alpha, N) returns the aerofoil surface
%   coordinates plus one far-field wake point for the Kutta condition.

NACAcode = char(NACAcode);
if numel(NACAcode) ~= 4 || ~all(isstrprop(NACAcode, 'digit'))
    error('NACAcode must be a four digit string, for example "2412".');
end

if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha)
    error('alpha must be a finite numeric scalar.');
end

if ~isscalar(N) || ~isnumeric(N) || ~isfinite(N) || N <= 0 || N ~= floor(N)
    error('N must be a positive integer.');
end

m = str2double(NACAcode(1)) / 100;
p = str2double(NACAcode(2)) / 10;
t = str2double(NACAcode(3:4)) / 100;

x = 1 - 0.5 * (1 - cos(2 * pi * (0:N) / N));
yc = zeros(size(x));
dydx = zeros(size(x));

if m > 0 && p > 0
    forward = x < p;
    aft = ~forward;

    yc(forward) = m / p^2 * (2 * p * x(forward) - x(forward).^2);
    dydx(forward) = 2 * m / p^2 * (p - x(forward));

    yc(aft) = m / (1 - p)^2 * (1 - 2 * p + 2 * p * x(aft) - x(aft).^2);
    dydx(aft) = 2 * m / (1 - p)^2 * (p - x(aft));
end

yt = 5 * t * ( ...
    0.2969 * sqrt(x) ...
    - 0.1260 * x ...
    - 0.3516 * x.^2 ...
    + 0.2843 * x.^3 ...
    - 0.1015 * x.^4);
theta = atan(dydx);

xUpper = x - yt .* sin(theta);
xLower = x + yt .* sin(theta);
zUpper = yc + yt .* cos(theta);
zLower = yc - yt .* cos(theta);

leadingEdge = find(xLower == min(xLower), 1);
X = [xLower(1:leadingEdge), xUpper(leadingEdge + 1:end), 1];
Z = [zLower(1:leadingEdge), zUpper(leadingEdge + 1:end), 0];

% Wake panel.
X(N + 2) = 1e8;
Z(N + 2) = X(end) * tand(alpha);
end
