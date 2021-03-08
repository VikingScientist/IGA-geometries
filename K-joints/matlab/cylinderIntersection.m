function [X ev_v] = cylinderIntersection(R, r, n, phi, theta, range),
% function [X ev_v] = cylinderIntersection(R, r, n, phi, theta, range),
% function [x y z] = cylinderIntersection(R, r, n, phi, theta),
%   returns sample points of the intersection curve between two cylinders. One
%   major cylinder with larger radius which is along the x-axis, while a smaller
%   cylinder with an angle of attack coming throught the large one. phi=theta=0
%   corresponds to the small cylinder along the y-axik.
%
%   parameters:
%      R      - radius of the large cylinder
%      r      - radius of the small cylinder (r<R)
%      n      - number of evaluation points
%      phi    - angle of attack in the xz-plane (along the major pipe). Range [0,pi/2)
%               default 0
%      theta  - angle of attack in the yz-plane (across the major pipe). Range: [0,pi/2)
%               default 0
%      range  - [start stop] angle sweep as seen from the smallest cylinder. Note that
%               start=0 refers to the points (x,y) = (1,0), and contiuing clockwise
%               default [pi/2, 5/2*pi]
%   returns
%      X       - nx3 matrix of coordinates for n evaluation points
%      ev_v    - evaluation angles (on the major cylinder)
% 

if(nargin<4)
	phi = 0;
end
if(nargin<5)
	theta = 0;
end
if(nargin<6)
	range = [pi/2, 5/2*pi];
end
if(R<r)
	disp 'Error: Must have R>r';
	return;
end
if theta<0 || theta >= pi/2,
	disp 'Error: Theta in range [0,pi/2)';
	return;
end
if phi<0 || phi >= pi/2,
	disp 'Error: Phi in range [0,pi/2)';
	return
end

h = (range(2)-range(1)) / (n-1);
u = [range(1):h:range(2)]';
% x = r.*sin(phi).*cos(u);
% y = r.*sin(u);
% z = sqrt(R.^2 - r.^2.*(sin(u).^2));

% declaring the rotation matrices
Rxz = eye(3);
Ryz = eye(3);
Rxz([1,3],[1,3]) = [cos(phi)  , -sin(phi);
                    sin(phi)  ,  cos(phi)];
Ryz([2,3],[2,3]) = [cos(theta), -sin(theta);
                    sin(theta),  cos(theta)];

% these are taken from maple worksheet (should lie in the same folder as this file)
x = cos(phi) .* r .* cos(u) - sin(phi) .* (-sin(phi) .* r .* cos(u) + sqrt(sin(phi) .^ 2 .* r .^ 2 .* cos(u) .^ 2 - r .^ 2 + R .^ 2 + r .^ 2 .* cos(u) .^ 2 .* cos(phi) .^ 2)) / cos(phi);
y = -sin(theta) .* sin(phi) .* r .* cos(u) + cos(theta) .* r .* sin(u) - sin(theta) .* (-sin(phi) .* r .* cos(u) + sqrt(sin(phi) .^ 2 .* r .^ 2 .* cos(u) .^ 2 - r .^ 2 + R .^ 2 + r .^ 2 .* cos(u) .^ 2 .* cos(phi) .^ 2));
z = R .* sqrt(0.1e1 - (sin(theta) .* sin(phi) .* r .* cos(u) - cos(theta) .* r .* sin(u) + sin(theta) .* (-sin(phi) .* r .* cos(u) + sqrt(sin(phi) .^ 2 .* r .^ 2 .* cos(u) .^ 2 - r .^ 2 + R .^ 2 + r .^ 2 .* cos(u) .^ 2 .* cos(phi) .^ 2))) .^ 2 / R .^ 2);

%%%% % plotting commands for debugging
% plot3(x,y,z, 'r-');
% xlabel 'x';
% ylabel 'y';
% zlabel 'z';

% hold on;
% plot3(-   R*ones(size(u)), R*cos(u), R*sin(u), 'b-');
% plot3(-.5*R*ones(size(u)), R*cos(u), R*sin(u), 'b-');
% plot3(  0*R*ones(size(u)), R*cos(u), R*sin(u), 'b-');
% plot3( .5*R*ones(size(u)), R*cos(u), R*sin(u), 'b-');
% plot3(    R*ones(size(u)), R*cos(u), R*sin(u), 'b-');

% Rot = Ryz*Rxz;
% plotX = Rot* [  r*cos(u), r*sin(u),  .7*R*ones(size(u))]';
% plot3(plotX(1,:), plotX(2,:) , plotX(3,:), 'g-');
% plotX = Rot* [  r*cos(u), r*sin(u),     R*ones(size(u))]';
% plot3(plotX(1,:), plotX(2,:) , plotX(3,:), 'g-');
% plotX = Rot* [  r*cos(u), r*sin(u), 1.5*R*ones(size(u))]';
% plot3(plotX(1,:), plotX(2,:) , plotX(3,:), 'g-');
% plotX = Rot* [  r*cos(u), r*sin(u), 2.1*R*ones(size(u))]';
% plot3(plotX(1,:), plotX(2,:) , plotX(3,:), 'g-');
% plotX = Rot* [  r*cos(u), r*sin(u), 3.3*R*ones(size(u))]';
% plot3(plotX(1,:), plotX(2,:) , plotX(3,:), 'g-');

% axis equal;
% hold off;
v  = pi - acos((sin(theta) .* sin(phi) .* r .* cos(u) - cos(theta) .* r .* sin(u) + sin(theta) .* (-sin(phi) .* r .* cos(u) + sqrt(sin(phi).^ 2 .* r.^ 2 .* cos(u).^ 2 - r.^ 2 + R.^ 2 + r.^ 2 .* cos(u).^ 2 .* cos(phi).^ 2))) / R);

X = [x, y, z];
ev_v = v;
