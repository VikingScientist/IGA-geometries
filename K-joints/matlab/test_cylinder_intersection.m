
r = 1;
R = 2;
p = 2;
theta = sin(r/R);
knot = [0,0,0, 1,1,1];

B = [r, 0, R;
     r, r, R;
	 0, r, sqrt(R^2-r^2)];
W = [1,          1,          1;
     sqrt(2)/2, sqrt(2)/2, sin(theta/2);
     1,          1,          1];

xi = 0:.05:1;
N  = getBSplineBasisAndDerivative(p, xi, knot);

W_tot = N'*W;

c = ( N'*(B.*W) ) ./ W_tot;

c(:,1).^2 + c(:,2).^2;
c(:,3).^2 + c(:,2).^2;

x = c(:,1);
y = c(:,2);
z = c(:,3);
plot3(x,y,z, 'r-');
axis([-.2 2.2 -.2 2.2 -.2 2.2]);
hold on;
plot3(B(:,1), B(:,2), B(:,3), 'bs ');
hold off;
close all;

% parameterizing the interesection curve (not B-splines), interpolating the points afterwards
n = 50;
h = 2*pi/ (n-1);
u = [0:h:2*pi]';
x = r*cos(u);
y = r*sin(u);
z = sqrt(R^2 - r^2*(sin(u).^2));

x.^2 + y.^2
z.^2 + y.^2

plot3(x,y,z, 'r-');
axis([-.2 2.2 -.2 2.2 -.2 2.2]);
