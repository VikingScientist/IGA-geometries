
L     =  40; % length of the large pipe in each direction  [-L, L]
l     =  40; % length of the small pipe in ONE direction   [0, l]
R     =  10; % radius of large pipe
r     =   7; % radius of small pipe
T     =   2; % thickness of large pipe
t     =   2; % thickness of small pipe
phi   =   pi/4; % angle of attack along the major pipe
theta =   pi/6; % angle of attack across the major pipe

center = [-R*atan(phi), R*atan(theta), R];

n = 20;                                      % number of interpolation points for each intersection HALF
p = 3;                                       % polynomial degree
knot = [zeros(1,p), 0:n-p, (n-p)*ones(1,p)]; % approximation space which we are looking for curves in
xi = knot(1):(knot(end)-knot(1))/(n-1):knot(end);

N = getBSplineBasisAndDerivative(p, xi, knot);

% we have four intersection curves since the different thickness pipes intersect at different curves
% for each of these curves, we split it into two pieces
X = cylinderIntersection(R-T/2, r-t/2, n, phi, theta, [pi/2, 3*pi/2]);
B_Inner_inner_left  = [N' \ X, ones(n,1)];
X = cylinderIntersection(R-T/2, r-t/2, n, phi, theta, [3*pi/2, 5*pi/2]);
B_Inner_inner_right = [N' \ X, ones(n,1)];

X = cylinderIntersection(R-T/2, r+t/2, n, phi, theta, [pi/2, 3*pi/2]);
B_Inner_outer_left  = [N' \ X, ones(n,1)];
X = cylinderIntersection(R-T/2, r+t/2, n, phi, theta, [3*pi/2, 5*pi/2]);
B_Inner_outer_right = [N' \ X, ones(n,1)];

X = cylinderIntersection(R+T/2, r-t/2, n, phi, theta, [pi/2, 3*pi/2]);
B_Outer_inner_left  = [N' \ X, ones(n,1)];
X = cylinderIntersection(R+T/2, r-t/2, n, phi, theta, [3*pi/2, 5*pi/2]);
B_Outer_inner_right = [N' \ X, ones(n,1)];

[X v] = cylinderIntersection(R+T/2, r+t/2, n, phi, theta, [pi/2, 3*pi/2]);
B_Outer_outer_left  = [N' \ X, ones(n,1)];
[X v] = cylinderIntersection(R+T/2, r+t/2, n, phi, theta, [3*pi/2, 5*pi/2]);
B_Outer_outer_right = [N' \ X, ones(n,1)];

% create the large pipe splitting curve
[X v] = majorCylinderSplit(R+T/2, r+t/2, n, phi, theta);
B_outer_split = [N' \ X, ones(n,1)];
[X v] = majorCylinderSplit(R-T/2, r+t/2, n, phi, theta);
B_inner_split = [N' \ X, ones(n,1)];

% Note thtat the sample angles v defined for the big pipe is different depending on which of the 4 intersection curves you're looking at

% create the (upper) end-curves for the large pipe
[X v] = cylinderIntersection(R-T/2, r+t/2, n, phi, theta, [pi/2, 3*pi/2]); % this is just to get the sample angles v
X = [(center(1)-L)*ones(size(v)), (R-T/2)*cos(v), (R-T/2)*sin(v)]; 
B_inner_upper_left_end = [N' \ X, ones(n,1)];
v = v(end:-1:1); % reversed direction on the right side
X = [(center(1)+L)*ones(size(v)), (R-T/2)*cos(v), (R-T/2)*sin(v)];
B_inner_upper_right_end = [N' \ X, ones(n,1)];
[X v] = cylinderIntersection(R+T/2, r+t/2, n, phi, theta, [pi/2, 3*pi/2]); % this is just to get the sample angles v
X = [(center(1)-L)*ones(size(v)), (R+T/2)*cos(v), (R+T/2)*sin(v)];
B_outer_upper_left_end = [N' \ X, ones(n,1)];
v = v(end:-1:1); % reversed direction on the right side
X = [(center(1)+L)*ones(size(v)), (R+T/2)*cos(v), (R+T/2)*sin(v)];
B_outer_upper_right_end = [N' \ X, ones(n,1)];

% using the same major-pipe-sweep-angles v to create the bottom end curves
[X v] = majorCylinderSplit(R+T/2, r+t/2, n, phi, theta); % this is just to get the sample angles v
X = [(center(1)-L)*ones(size(v)), (R+T/2)*cos(v), (R+T/2)*sin(v)];
B_outer_lower_left_end = [N' \ X, ones(n,1)];
X = [(center(1)+L)*ones(size(v)), (R+T/2)*cos(v), (R+T/2)*sin(v)];
B_outer_lower_right_end = [N' \ X, ones(n,1)];
[X v] = majorCylinderSplit(R-T/2, r+t/2, n, phi, theta);% this is just to get the sample angles v
X = [(center(1)-L)*ones(size(v)), (R-T/2)*cos(v), (R-T/2)*sin(v)];
B_inner_lower_left_end = [N' \ X, ones(n,1)];
X = [(center(1)+L)*ones(size(v)), (R-T/2)*cos(v), (R-T/2)*sin(v)];
B_inner_lower_right_end = [N' \ X, ones(n,1)];

% create the small pipe end points
Rxz = eye(3);
Ryz = eye(3);
Rxz([1,3],[1,3]) = [cos(phi)  , -sin(phi);
                    sin(phi)  ,  cos(phi)];
Ryz([2,3],[2,3]) = [cos(theta), -sin(theta);
                    sin(theta),  cos(theta)];
Rot = Ryz*Rxz;
v = [pi/2:pi/(n-1):3*pi/2]';
X = [(r+t/2)*cos(v), (r+t/2)*sin(v), l*ones(size(v))] + [ones(n,1)*center];
X = [Rot * X']';
B_outer_small_left_end = [N' \ X, ones(n,1)];
X = [(r-t/2)*cos(v), (r-t/2)*sin(v), l*ones(size(v))] + [ones(n,1)*center];
X = [Rot * X']';
B_inner_small_left_end = [N' \ X, ones(n,1)];
v = v + pi;
X = [(r+t/2)*cos(v), (r+t/2)*sin(v), l*ones(size(v))] + [ones(n,1)*center];
X = [Rot * X']';
B_outer_small_right_end = [N' \ X, ones(n,1)];
X = [(r-t/2)*cos(v), (r-t/2)*sin(v), l*ones(size(v))] + [ones(n,1)*center];
X = [Rot * X']';
B_inner_small_right_end = [N' \ X, ones(n,1)];


curves = cell(14,1);
curves{1}  = B_Inner_inner_left;
curves{2}  = B_Inner_inner_right;
curves{3}  = B_Inner_outer_left;
curves{4}  = B_Inner_outer_right;
curves{5}  = B_Outer_inner_left;
curves{6}  = B_Outer_inner_right;
curves{7}  = B_Outer_outer_left;
curves{8}  = B_Outer_outer_right;
curves{9}  = B_outer_split;
curves{10} = B_inner_split;
curves{11} = B_outer_lower_left_end;
curves{12} = B_inner_lower_left_end;
curves{13} = B_outer_lower_right_end;
curves{14} = B_inner_lower_right_end;
curves{15} = B_outer_upper_left_end;
curves{16} = B_inner_upper_left_end;
curves{17} = B_outer_upper_right_end;
curves{18} = B_inner_upper_right_end;
curves{19} = B_outer_small_right_end;
curves{20} = B_inner_small_right_end;
curves{21} = B_outer_small_left_end;
curves{22} = B_inner_small_left_end;

plotN = 100;
xi = knot(1):(knot(end)-knot(1))/(plotN-1):knot(end);
N = getBSplineBasisAndDerivative(p, xi, knot);
figure; hold on;
for i=1:length(curves),
	C = N'*curves{i};
	W = C(:,4);
	C(:,1) = C(:,1) ./ W;
	C(:,2) = C(:,2) ./ W;
	C(:,3) = C(:,3) ./ W;
	plot3(C(:,1), C(:,2), C(:,3), 'b-');
end
axis equal;
hold off;

disp 'press any key to create surfaces';
pause;

p = [p 1];
knot = struct('xi', knot, 'eta', [0,0,1,1]);
B_surf_outer_left_top = zeros(n, 2, 3);
B_surf_outer_left_top(:,1,:) = B_outer_upper_left_end(:,1:3);
B_surf_outer_left_top(:,2,:) = B_Outer_outer_left(:,1:3);
B_surf_outer_right_top = zeros(n, 2, 3);
B_surf_outer_right_top(:,1,:) = B_outer_upper_right_end(:,1:3);
B_surf_outer_right_top(:,2,:) = B_Outer_outer_right(:,1:3);

B_surf_inner_left_top = zeros(n, 2, 3);
B_surf_inner_left_top(:,1,:) = B_inner_upper_left_end(:,1:3);
B_surf_inner_left_top(:,2,:) = B_Inner_outer_left(:,1:3);
B_surf_inner_right_top = zeros(n, 2, 3);
B_surf_inner_right_top(:,1,:) = B_inner_upper_right_end(:,1:3);
B_surf_inner_right_top(:,2,:) = B_Inner_outer_right(:,1:3);

B_surf_outer_left_bottom = zeros(n, 2, 3);
B_surf_outer_left_bottom(:,1,:) = B_outer_lower_left_end(:,1:3);
B_surf_outer_left_bottom(:,2,:) = B_outer_split(:,1:3);
B_surf_outer_right_bottom = zeros(n, 2, 3);
B_surf_outer_right_bottom(:,1,:) = B_outer_lower_right_end(:,1:3);
B_surf_outer_right_bottom(:,2,:) = B_outer_split(:,1:3);

B_surf_inner_left_bottom = zeros(n, 2, 3);
B_surf_inner_left_bottom(:,1,:) = B_inner_lower_left_end(:,1:3);
B_surf_inner_left_bottom(:,2,:) = B_inner_split(:,1:3);
B_surf_inner_right_bottom = zeros(n, 2, 3);
B_surf_inner_right_bottom(:,1,:) = B_inner_lower_right_end(:,1:3);
B_surf_inner_right_bottom(:,2,:) = B_inner_split(:,1:3);

B_surf_inner_left_intersection = zeros(n, 2, 3);
B_surf_inner_left_intersection(:,1,:) = B_Inner_outer_left(:,1:3);
B_surf_inner_left_intersection(:,2,:) = B_Inner_inner_left(:,1:3);
B_surf_inner_right_intersection = zeros(n, 2, 3);
B_surf_inner_right_intersection(:,1,:) = B_Inner_outer_right(:,1:3);
B_surf_inner_right_intersection(:,2,:) = B_Inner_inner_right(:,1:3);
B_surf_outer_left_intersection = zeros(n, 2, 3);
B_surf_outer_left_intersection(:,1,:) = B_Outer_outer_left(:,1:3);
B_surf_outer_left_intersection(:,2,:) = B_Outer_inner_left(:,1:3);
B_surf_outer_right_intersection = zeros(n, 2, 3);
B_surf_outer_right_intersection(:,1,:) = B_Outer_outer_right(:,1:3);
B_surf_outer_right_intersection(:,2,:) = B_Outer_inner_right(:,1:3);

B_surf_small_outer_left = zeros(n, 2, 3);
B_surf_small_outer_left(:,1,:) = B_outer_small_left_end(:,1:3);
B_surf_small_outer_left(:,2,:) = B_Outer_outer_left(:,1:3);
B_surf_small_inner_left = zeros(n, 2, 3);
B_surf_small_inner_left(:,1,:) = B_inner_small_left_end(:,1:3);
B_surf_small_inner_left(:,2,:) = B_Outer_inner_left(:,1:3);
B_surf_small_inner_right = zeros(n, 2, 3);
B_surf_small_inner_right(:,1,:) = B_inner_small_right_end(:,1:3);
B_surf_small_inner_right(:,2,:) = B_Outer_inner_right(:,1:3);
B_surf_small_outer_right = zeros(n, 2, 3);
B_surf_small_outer_right(:,1,:) = B_outer_small_right_end(:,1:3);
B_surf_small_outer_right(:,2,:) = B_Outer_outer_right(:,1:3);

% shut the hole tight by filling in the small pipe
B_surf_fill_top = zeros(n, 2, 3);
B_surf_fill_top(:,1,:) = B_Outer_inner_left(:,1:3);
B_surf_fill_top(:,2,:) = B_Outer_inner_right(end:-1:1,1:3);
B_surf_fill_bottom = zeros(n, 2, 3);
B_surf_fill_bottom(:,1,:) = B_Inner_inner_left(:,1:3);
B_surf_fill_bottom(:,2,:) = B_Inner_inner_right(end:-1:1,1:3);

surfaces = cell(1);

surfaces{1}  = B_surf_outer_left_top;
surfaces{2}  = B_surf_outer_right_top;
surfaces{3}  = B_surf_inner_left_top;
surfaces{4}  = B_surf_inner_right_top;
surfaces{5}  = B_surf_outer_left_bottom;
surfaces{6}  = B_surf_outer_right_bottom;
surfaces{7}  = B_surf_inner_left_bottom;
surfaces{8}  = B_surf_inner_right_bottom;
surfaces{9}  = B_surf_inner_right_intersection;
surfaces{10} = B_surf_inner_left_intersection;
surfaces{11} = B_surf_outer_right_intersection;
surfaces{12} = B_surf_outer_left_intersection;
surfaces{13} = B_surf_small_outer_left;
surfaces{14} = B_surf_small_inner_left;
surfaces{15} = B_surf_small_outer_right;
surfaces{16} = B_surf_small_inner_right;
surfaces{17} = B_surf_fill_top;
surfaces{18} = B_surf_fill_bottom;

plotX = 40;
plotY = 3;
Nx = getBSplineBasisAndDerivative(p(1), knot.xi(1):(knot.xi(end)-knot.xi(1))/(plotX-1):knot.xi(end), knot.xi);
Ny = getBSplineBasisAndDerivative(p(2), 0:1/(plotY-1):1, knot.eta);
figure; hold on;
for i=1:length(surfaces),
	X = Nx'*surfaces{i}(:,:,1)*Ny;
	Y = Nx'*surfaces{i}(:,:,2)*Ny;
	Z = Nx'*surfaces{i}(:,:,3)*Ny;
	surf(X, Y, Z);
end
axis equal;
hold off;

B_vol_left_bottom = zeros(n,2,2,3);
B_vol_left_bottom(:,:,1,:) = B_surf_inner_left_bottom;
B_vol_left_bottom(:,:,2,:) = B_surf_outer_left_bottom;
B_vol_left_top = zeros(n,2,2,3);
B_vol_left_top(:,:,1,:) = B_surf_inner_left_top;
B_vol_left_top(:,:,2,:) = B_surf_outer_left_top;
B_vol_right_bottom = zeros(n,2,2,3);
B_vol_right_bottom(:,:,1,:) = B_surf_inner_right_bottom;
B_vol_right_bottom(:,:,2,:) = B_surf_outer_right_bottom;
B_vol_right_top = zeros(n,2,2,3);
B_vol_right_top(:,:,1,:) = B_surf_inner_right_top;
B_vol_right_top(:,:,2,:) = B_surf_outer_right_top;
B_vol_left_intersection = zeros(n,2,2,3);
B_vol_left_intersection(:,:,1,:) = B_surf_inner_left_intersection;
B_vol_left_intersection(:,:,2,:) = B_surf_outer_left_intersection;
B_vol_right_intersection = zeros(n,2,2,3);
B_vol_right_intersection(:,:,1,:) = B_surf_inner_right_intersection;
B_vol_right_intersection(:,:,2,:) = B_surf_outer_right_intersection;
B_vol_small_right = zeros(n,2,2,3);
B_vol_small_right(:,:,1,:) = B_surf_small_inner_right;
B_vol_small_right(:,:,2,:) = B_surf_small_outer_right;
B_vol_small_left = zeros(n,2,2,3);
B_vol_small_left(:,:,1,:) = B_surf_small_inner_left;
B_vol_small_left(:,:,2,:) = B_surf_small_outer_left;
B_vol_fill = zeros(n,2,2,3);
B_vol_fill(:,:,1,:) = B_surf_fill_top;
B_vol_fill(:,:,2,:) = B_surf_fill_bottom;

volumes = cell(1);
volumes{1} = B_vol_left_bottom;
volumes{2} = B_vol_left_top;
volumes{3} = B_vol_right_bottom;
volumes{4} = B_vol_right_top;
volumes{5} = B_vol_left_intersection;
volumes{6} = B_vol_right_intersection;
volumes{7} = B_vol_small_left;
volumes{8} = B_vol_small_right;
volumes{9} = B_vol_fill;

knot = struct('xi', knot.xi, 'eta', knot.eta, 'zeta', [0,0,1,1]);
str = write_g2_volume(volumes{1}, knot);
for i=2:length(volumes),
	str = sprintf('%s\n\n%s', str, write_g2_volume(volumes{i}, knot));
end

fid = fopen('out.g2', 'w');
fprintf(fid, '%s', str);
fclose(fid);
