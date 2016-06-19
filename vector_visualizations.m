function [] = vector_visualizations(R_2_to_W, R_W_to_1,cam1,cam2,n,p)% sample data from first point in our c++ run
R1 = R_W_to_1;
R2 = transpose(R_2_to_W);
normal = n;
center1 = cam1;
center2 = cam2;
point = p;


%% here we go, lets start plotting
figure(1)
scale = 2; % scaling how big the vector arrows should look like
xyz_world = [0;0;0]; % world frame origin

%world coordinate frame
% qworldx = quiver3(0,0,0,1,0,0);
hold on
% qworldx.Color = 'k';
% qworldx.AutoScaleFactor = scale;
% qworldy = quiver3(0,0,0,0,1,0);
% qworldy.Color = 'k';
% qworldy.AutoScaleFactor = scale;
% qworldz = quiver3(0,0,0,0,0,1);
% qworldz.Color = 'r';
% qworldz.AutoScaleFactor = scale;

% camera 1 coordinate frame
qcamera1x = quiver3(center1(1),center1(2),center1(3),R1(1,1),R1(1,2),R1(1,3));
qcamera1x.Color = 'b';
qcamera1x.AutoScaleFactor = scale;
qcamera1y = quiver3(center1(1),center1(2),center1(3),R1(2,1),R1(2,2),R1(2,3));
qcamera1y.Color = 'b';
qcamera1y.AutoScaleFactor = scale;
qcamera1z = quiver3(center1(1),center1(2),center1(3),R1(3,1),R1(3,2),R1(3,3));
qcamera1z.Color = 'g';
qcamera1z.AutoScaleFactor = scale;

% normal from the surface point
qnormal = quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3));
qnormal.Color = 'm';
qnormal.AutoScaleFactor = scale;
qnormal.LineWidth = 3; %make it fat!
pnormal = scatter3(point(1),point(2),point(3),30,'r');

% camera 2 coordinate frame
qcamera2x = quiver3(center2(1),center2(2),center2(3),R2(1,1),R2(1,2),R2(1,3));
qcamera2x.Color = 'b';
qcamera2x.AutoScaleFactor = scale;
qcamera2y = quiver3(center2(1),center2(2),center2(3),R2(2,1),R2(2,2),R2(2,3));
qcamera2y.Color = 'b';
qcamera2y.AutoScaleFactor = scale;
qcamera2z = quiver3(center2(1),center2(2),center2(3),R2(3,1),R2(3,2),R2(3,3));
qcamera2z.Color = 'c';
qcamera2z.AutoScaleFactor = scale;

% equal unit size on all axis
axis equal
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');

pause
clf