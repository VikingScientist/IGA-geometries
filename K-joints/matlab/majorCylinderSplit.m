function [X ev_v] = majorCylinderSplit(R, r, n, phi, theta, range),
% function [X ev_v] = majorCylinderSplit(R, r, n, phi, theta, range),
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
	range = [pi/2, 3/2*pi];
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

% this is taken from maple - its the start and stop angle (at the intersection points) as
% seen from the major cylinder
angle_back  = pi - acos((sin(theta) * sin(phi) * r * cos(range(1)) - cos(theta) * r * sin(range(1)) + sin(theta) * (-sin(phi) * r * cos(range(1)) + sqrt(sin(phi) ^ 2 * r ^ 2 * cos(range(1)) ^ 2 - r ^ 2 + R ^ 2 + r ^ 2 * cos(range(1)) ^ 2 * cos(phi) ^ 2))) / R);
angle_front = pi - acos((sin(theta) * sin(phi) * r * cos(range(2)) - cos(theta) * r * sin(range(2)) + sin(theta) * (-sin(phi) * r * cos(range(2)) + sqrt(sin(phi) ^ 2 * r ^ 2 * cos(range(2)) ^ 2 - r ^ 2 + R ^ 2 + r ^ 2 * cos(range(2)) ^ 2 * cos(phi) ^ 2))) / R);

% this is taken from maple - its the start and stop x-coordinate (at the intersection points)
x_back  = cos(phi) * r * cos(pi/2) - sin(phi) * (-sin(phi) * r * cos(pi/2) + sqrt(sin(phi) ^ 2 * r ^ 2 * cos(pi/2) ^ 2 - r ^ 2 + R ^ 2 + r ^ 2 * cos(pi/2) ^ 2 * cos(phi) ^ 2)) / cos(phi);
x_front = cos(phi) * r * cos(3*pi/2) - sin(phi) * (-sin(phi) * r * cos(3*pi/2) + sqrt(sin(phi) ^ 2 * r ^ 2 * cos(3*pi/2) ^ 2 - r ^ 2 + R ^ 2 + r ^ 2 * cos(3*pi/2) ^ 2 * cos(phi) ^ 2)) / cos(phi);

v = [angle_front:(2*pi+angle_back-angle_front)/(n-1):(2*pi+angle_back)]';

x = [x_front:(x_back-x_front)/(n-1):x_back]';
y = R*cos(v);
z = R*sin(v);

X = [x, y, z];
ev_v = v;
