function str = write_g2_volume(B, knot),

[n1 n2 n3 dim] = size(B);

if dim ~=3,
	disp 'Error: B must be of size n1 x n2 x n3 x 3'
	return;
end

str = sprintf('700 1 0 0\n3 0');
str = strcat(str, sprintf('\n%d %d', n1, length(knot.xi)-n1));
knotxi = sprintf('%.16f ', knot.xi);
str = strcat(str, sprintf('\n%s', knotxi));
str = strcat(str, sprintf('\n%d %d', n2, length(knot.eta)-n2));
knotxi = sprintf('%.16f ', knot.eta);
str = strcat(str, sprintf('\n%s', knotxi));
str = strcat(str, sprintf('\n%d %d', n3, length(knot.zeta)-n3));
knotxi = sprintf('%.16f ', knot.zeta);
str = strcat(str, sprintf('\n%s', knotxi));
for k=1:n3
	for j=1:n2,
		for i=1:n1,
			oneComp = sprintf('%.16f ', B(i,j,k,:));
			str = strcat(str, sprintf('\n%s', oneComp));
		end
	end
end
str = strcat(str, sprintf('\n'));
