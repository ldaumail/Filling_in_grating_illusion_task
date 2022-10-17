function [x, y] = drawEllipse(w,h,xcenter,ycenter)
%w = width
%h = height
%w = figureDiameterhorz;
%h =figureDiameterVert;
%Loic Daumail using code from Roger Stafford and Image Analyst, 10/20/2021
% Define parameters.
x1 = xcenter;
x2 = xcenter;
y1 = ycenter-h/2;
y2 = ycenter+h/2;
eccentricity = w/h;
numPoints = 300; % Less for a coarser ellipse, more for a finer resolution.
% Make equations:
a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
b = a * eccentricity; %sqrt(1-eccentricity^2);
t = linspace(0, 2 * pi, numPoints); % Absolute angle parameter
X = a * cos(t);
Y = b * sin(t);
% Compute angles relative to (x1, y1).
angles = atan2(y2 - y1, x2 - x1);
x = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
y = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
% Plot the ellipse as a blue curve.
%plot(x,y,'b-', 'LineWidth', 2);	% Plot ellipse
%grid on;
%axis equal

end