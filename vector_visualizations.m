function [] = vector_visualizations(R_2_to_W, R_W_to_1,cam1,cam2,n,p)% sample data from first point in our c++ run

% % sample data from first point in our c++ run
% R_W_to_1 = [0.7578819431221295, -0.4160528855759703, -0.5025089488081022;
%  0.283021121067612, 0.9036777088471941, -0.3213498991998712;
%  0.5878046378329685, 0.1013247703675849, 0.8026326884421708];
% 
% R_W_to_2 = [-0.03329457833514037, 0.9994455818369928, 0;
%  -0.7337150549710932, -0.02444228461994558, 0.679017520268309;
%  0.678641, 0.0226076, 0.7341220000000001];
% 
% cam1 = [0.8181166517939999;
%  2.77555756156e-17;
%  0.15391315114];
% 
% cam2 = [2.770106104830846;
%  1.149446642778977;
%  3.806435863614157];
% 
% n = [-0.678641;
%  -0.0226076;
%  -0.7341220000000001];
% 
% p = [2.77106;1.14948;3.80746]; 

R_W_to_1 = [0.9853917887504152, -0.1520612684086295, -0.07668366306956452;
 0.1428960436721329, 0.9832135069660647, -0.1134544398834716;
 0.09264844431861587, 0.1008392774487159, 0.9905794768721412];

R_W_to_2 = [-0.9120278855581706, -0.410128194549329, 0;
 0.3660446254729072, -0.8139964777521151, -0.4510222459807821;
 0.184977, -0.411345, 0.892513];

n = [-0.184977;
 0.411345;
 -0.892513];

cam1 = [0.7912039016370001;
 -0.0115708260483;
 0.137853001671];

cam2 = [2.672833584137161;
 1.098753351287903;
 3.83942062696008];

p = [2.67306;1.09825;3.8405];


R1 = R_W_to_1;
R2 = R_W_to_2; %transpose(R_2_to_W);
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
grid on
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