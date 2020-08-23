figure(1);
fileID = fopen("initialized.txt");
C = textscan(fileID, "%f %f %f %f %f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
u1 = cell2mat(C(1,3));
u2 = -cell2mat(C(1,4));
v1 = cell2mat(C(1,5));
v2 = -cell2mat(C(1,6));
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(2)
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");

figure(3)
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
fileID = fopen("discontinuities.txt");
C = textscan(fileID, "%f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
scatter(xs,ys,"filled");


figure(4);
fileID = fopen("z1z2.txt");
C = textscan(fileID, "%f %f %f %f %f %f");
xs = cell2mat(C(1,1));
ys = cell2mat(C(1,2));
u1 = cell2mat(C(1,3));
u2 = cell2mat(C(1,4));
v1 = cell2mat(C(1,5));
v2 = cell2mat(C(1,6));
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(5)
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
%quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(7);
fileID = fopen("ab.txt");
C = textscan(fileID, "%f %f %f %f %f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
u1 = cell2mat(C(1,3));
u2 = -cell2mat(C(1,4));
v1 = cell2mat(C(1,5));
v2 = -cell2mat(C(1,6));
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
%quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
%quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(8)
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
%quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");


figure(4);
fileID = fopen("optimized.txt");
C = textscan(fileID, "%f %f %f %f %f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
u1 = cell2mat(C(1,3));
u2 = -cell2mat(C(1,4));
v1 = cell2mat(C(1,5));
v2 = -cell2mat(C(1,6));
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(5)
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(6)
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");

figure(9);
fileID = fopen("edge_vertices.txt");
C = textscan(fileID, "%f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
scatter(xs,ys);

figure(10);
fileID = fopen("constraints.txt");
C = textscan(fileID, "%f %f %f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
u1 = cell2mat(C(1,3));
u2 = -cell2mat(C(1,4));
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");

figure(11)
fileID = fopen("magnitudes.txt");
C = textscan(fileID, "%f %f %f %f %f %f %f %f");
xs = cell2mat(C(1,1));
ys = -cell2mat(C(1,2));
u1 = cell2mat(C(1,3));
u2 = -cell2mat(C(1,4));
v1 = cell2mat(C(1,5));
v2 = -cell2mat(C(1,6));
lambda_u = cell2mat(C(1,7));
lambda_v = cell2mat(C(1,8));
u_neg = lambda_u < 0;
u_pos = lambda_u >= 0;
v_neg = lambda_v < 0;
v_pos = lambda_v >= 0;
quiver(xs(u_neg), ys(u_neg), lambda_u(u_neg).*u1(u_neg), lambda_u(u_neg).*u2(u_neg), .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
quiver(xs(u_pos), ys(u_pos), lambda_u(u_pos).*u1(u_pos), lambda_u(u_pos).*u2(u_pos), .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
quiver(xs(v_neg), ys(v_neg), lambda_v(v_neg).*v1(v_neg), lambda_v(v_neg).*v2(v_neg), .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
quiver(xs(v_pos), ys(v_pos), lambda_v(v_pos).*v1(v_pos), lambda_v(v_pos).*v2(v_pos), .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");

