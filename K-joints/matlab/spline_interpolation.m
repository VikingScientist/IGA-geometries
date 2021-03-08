function [c knot N] = spline_interpolation(p, n, f, knot, xi),
% function [c knot] = spline_interpolation(p, n, f, knot, xi),
%     Computes a spline interpolation. Note that this is a 1D function interpolation,
%     and that the boundary conditions are set to 'free' (i.e. no restrictions on boundary).
%     Also generates a plot of the approximation compared to the exact function
%     parameters:
%         p      - the polynomial order of the basis
%         n      - the number of control points requested (ignored if knot is specified)
%         f      - function which is to be interpolated
%         [knot] - [optional] knot vector describing spline space (default uniform)
%         [xi]   - [optional] interpolation points (default unifrom including endpoints)
%     returns:
%         c      - control points which interpolates the requested function
%         knot   - knot vector describing spline space
%         N      - the matrix which is used for calculating the control points (N \ f(xi) )
%     example:
%         % computes a cubic spline interpolation of sinus on the interval [0,8*pi] using 20 pts
%         c = spline_interpolation(3, 20, @sin);
%
%         % computes an interpolation of the heavyside function
%         p = 4;  % degere
%         n = 10; % interpolation pts;
%         f = @(x) x<4*pi;
%         c = spline_interpolation(p,n,f);


if nargin < 3,
	disp 'Error - not enough input arguments';
	disp '        see documentation';
	return;
end
if nargin < 4,
	knot = zeros(1, n+p+1);
	knot(1:p) = 0;
	knot(p+1:n+1) = 0:8*pi/(n-p):8*pi;
	knot(n+2:end) = 8*pi;
end
if nargin < 5,
	n = length(knot) - p - 1; % knot vector giving number of control points, not input-n
	h = (knot(end)-knot(1))/(n-1);
	xi = knot(1):h:knot(end);
end

xi = xi(:); % makes xi a row-vector
N = getBSplineBasisAndDerivative(p, xi, knot);
N = N';
c = N \ f(xi);

plotN = 1000;
h = (knot(end)-knot(1))/(plotN+1);
x = knot(1):h:knot(end);
Np = getBSplineBasisAndDerivative(p, x, knot);
figure; plot(x, Np'*c, 'b-');
hold on;
plot(x, f(x), 'r-');
plot(xi, f(xi), 'bo ');
legend('estimated', 'exact', 'interp.pts');
hold off;

figure;
plot(x, Np'*c-f(x'));
title('error');
