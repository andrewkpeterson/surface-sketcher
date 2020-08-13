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
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");

figure(3)
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

%{
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

figure(8)
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "on");
hold on;
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "on");
%quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
%}

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
quiver(xs, ys, v1, v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -v1, -v2, .3, "LineWidth", 2, "Color", [0,0,1], "ShowArrowHead", "off");

figure(6)
quiver(xs, ys, u1, u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");
hold on;
quiver(xs, ys, -u1, -u2, .3, "LineWidth", 2, "Color", [1,0,0], "ShowArrowHead", "off");