% sample data from first point in our c++ run
R_W_to_1 = [0.7578819431221295, -0.4160528855759703, -0.5025089488081022;
 0.283021121067612, 0.9036777088471941, -0.3213498991998712;
 0.5878046378329685, 0.1013247703675849, 0.8026326884421708];

R_W_to_2 = [-0.03329457833514037, 0.9994455818369928, 0;
 -0.7337150549710932, -0.02444228461994558, 0.679017520268309;
 -0.678641, -0.0226076, -0.7341220000000001];

C1 = [0.8181166517939999;
 2.77555756156e-17;
 0.15391315114];

C2 = [2.770106104830846;
 1.149446642778977;
 3.806435863614157];

normal_W = [-0.678641;
 -0.0226076;
 -0.7341220000000001];

point_W = [2.77106;1.14948;3.80746]; 

%% parameter selection. Useless thing, just did it to be able to quickly
%switch between different matrixes from above.
R1 = R_W_to_1;
R2 = R_W_to_2;
normal = normal_W;
center1 = C1;
center2 = C2;
point = point_W;


%% here we go, lets start plotting
figure(1)
scale = 2; % scaling how big the vector arrows should look like
xyz_world = [0;0;0]; % world frame origin

%world coordinate frame
qworldx = quiver3(0,0,0,1,0,0);
hold on
qworldx.Color = 'k';
qworldx.AutoScaleFactor = scale;
qworldy = quiver3(0,0,0,0,1,0);
qworldy.Color = 'k';
qworldy.AutoScaleFactor = scale;
qworldz = quiver3(0,0,0,0,0,1);
qworldz.Color = 'r';
qworldz.AutoScaleFactor = scale;

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