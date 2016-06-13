%% Circle generation
sz=400;
rad=100;
clear RGB
RGB(1:sz,1:sz,1:3)=255;
[x y]= find(RGB==255);
xc=ceil((sz+1)/2);
yc=ceil((sz+1)/2);

r=rad.^2;
d=find(((x-xc).^2+(y-yc).^2) <= r);
for i=1:size(d,1)
      
    RGB(x(d(i)),y(d(i)),:)=0;
      
end

B=rgb2gray(RGB);

ED=edge(B);
SE=strel('disk',1);
cir=~imdilate(ED,SE);

cir(3:10,3:10) = 0;

figure(1),imshow(cir);
imwrite(cir,'circle_original.jpeg');



    %% homogenious transfor
clear image_final;
theta = 50;
% R = [cosd(theta) -sind(theta) 0;
%      sind(theta) cosd(theta)  0;
%      0          0           1];

R = [1 0 0;
    0 cosd(theta) -sind(theta);
    0 sind(theta) cosd(theta)];

T = [-200;-200;400];

K = [400 0 200;
     0 400 200;
     0 0 1];
H = K*[R,T];


for i=1:1:400
    for j=1:1:400
        z = 0;
        point = [i;j;z;1];
        point_new = H * point;
        u_coord = round(point_new(1)/point_new(3));
        v_coord = round(point_new(2)/point_new(3));
        %if((u_coord > 0 && u_coord < 400) && (v_coord > 0 && v_coord < 400))
        image_final(u_coord,v_coord) = cir(i,j);
        %end
        
    end
end
figure(2)
imshow(image_final);
imwrite(image_final,'circle_warped.jpeg');
            

%% circle camera view
test=0;
if test
clear image_final;
theta = -50;
% R = [cosd(theta) -sind(theta) 0;
%      sind(theta) cosd(theta)  0;
%      0          0           1];

R = [1 0 0;
    0 cosd(theta) -sind(theta);
    0 sind(theta) cosd(theta)];

T = [200;200;400];

K = [400 0 200;
     0 400 200;
     0 0 1];
H = K*[R,T];


for i=1:1:400
    for j=1:1:400
        z = 0;
        point = [i;j;z;1];
        point_new = H * point;
        u_coord = round(point_new(1)/point_new(3));
        v_coord = round(point_new(2)/point_new(3));
        %if((u_coord > 0 && u_coord < 400) && (v_coord > 0 && v_coord < 400))
        image_final(u_coord,v_coord) = cir(i,j);
        %end
        
    end
end
figure(3)
imshow(image_final);
%imwrite(image_final,'/Users/markopanjek/Documents/DOKUMENTI/Sola/ETH/3D_Vision/semester_project/circle_warped.jpeg');
end
            
%% UNWARPING USING HOMOGRAPHY FROM C++ (Torstens meeting, 8.6.2016)
clear image_unwarped;

H = [1, -0.373383055017823, 0.0085315077195105;
 0, 0.9748836410988672, -0.01861046434370905;
 0, -0.001866915275089115, 1.000042657538598];

[maxX,maxY]= size(image_final);
for i=1:1:400
    for j=1:1:400
        point = [i;j;1];
        point_new = inv(H) * point;
        u_coord = round(point_new(1)/point_new(3));
        v_coord = round(point_new(2)/point_new(3));
        image_unwarped(i,j) = image_final(u_coord,v_coord);
    end
end
figure(4)
imshow(image_unwarped);


%% WARPING A COLOR IMAGE

clear color_warped;
color_original = imread('color_circle.jpg');

theta = 50;
% R = [cosd(theta) -sind(theta) 0;
%      sind(theta) cosd(theta)  0;
%      0          0           1];

R = [1 0 0;
    0 cosd(theta) -sind(theta);
    0 sind(theta) cosd(theta)];

T = [-200;-200;400];

K = [400 0 200;
     0 400 200;
     0 0 1];
H = K*[R,T];

for i=1:1:400
    for j=1:1:400
        z = 0;
        point = [i;j;z;1];
        point_new = H * point;
        u_coord = round(point_new(1)/point_new(3));
        v_coord = round(point_new(2)/point_new(3));
        %if((u_coord > 0 && u_coord < 400) && (v_coord > 0 && v_coord < 400))
        color_warped(u_coord,v_coord,:) = color_original(i,j,:);
        %end
        
    end
end
figure(2)
imshow(color_warped);
imwrite(color_warped,'color_warped.jpeg');
 %%            
%lio
% unwarp of color_circle
clear image_unwarped;

H = [1, -0.373383055017823, 0.0085315077195105;
 0, 0.9748836410988672, -0.01861046434370905;
 0, -0.001866915275089115, 1.000042657538598];

for i=1:1:400
    for j=1:1:400
        point = [i;j;1];
        point_new = inv(H) * point;
        u_coord = round(point_new(1)/point_new(3));
        v_coord = round(point_new(2)/point_new(3));    
        image_unwarped(i,j,:) = color_warped(u_coord,v_coord,:);
    end
end
figure(3)
imshow(image_unwarped);

%% %% WARPING A COLOR IMAGE using bilinear interpolation

clear color_warped;
color_original = imread('color_circle.jpg');

theta = 50;
% R = [cosd(theta) -sind(theta) 0;
%      sind(theta) cosd(theta)  0;
%      0          0           1];

R = [1 0 0;
    0 cosd(theta) -sind(theta);
    0 sind(theta) cosd(theta)];

T = [-200;-200;400];

K = [400 0 200;
     0 400 200;
     0 0 1];
H = K*[R,T];

for u_coord= 1:1:400
    for v_coord=1:1:232
        z = 0;
        point = [u_coord;v_coord;z;1];
        point_new = H * point;
        i = round(point_new(1)/point_new(3));
        j = round(point_new(2)/point_new(3));
        if((i > 1 && i < 399) && (j > 1 && j < 399))
            
        x1 = i-1;
        x2 = i+1;
        y1 = j-1;
        y2 = j+1;
        x = i;
        y = j;
            
        f_Q11 = color_original(i-1, j-1,:);
        f_Q21 = color_original(i+1, j-1,:);
        f_Q12 = color_original(i-1, j+1,:);
        f_Q22 = color_original(i+1, j+1,:);
        
        f_xy1 = ((x2-x)/(x2-x1))* f_Q11 + ((x-x1)/(x2-x1))* f_Q21;
        f_xy2 = ((x2-x)/(x2-x1))* f_Q12 + ((x-x1)/(x2-x1))* f_Q22;
        
        f_xy = ((y2-y)/(y2-y1))* f_xy1 + ((y-y1)/(y2-y1))* f_xy2;
        
        color_warped(u_coord,v_coord,:) = f_xy; %color_original(i,j,:);
        end
        
    end
end
figure(2)
imshow(color_warped);
imwrite(color_warped,'color_warped.jpeg');
%%
%sonali

% unwarp of color_circle
clear image_unwarped;

H = [1, -0.373383055017823, 0.0085315077195105;
 0, 0.9748836410988672, -0.01861046434370905;
 0, -0.001866915275089115, 1.000042657538598];

for i=1:1:400
    for j=1:1:400
        point = [i;j;1];
        point_new = inv(H) * point;
        U_coord = point_new(1)/point_new(3);
        V_coord = point_new(2)/point_new(3);
        u_coord = floor(point_new(1)/point_new(3));
        v_coord = floor(point_new(2)/point_new(3));
      if((u_coord > 1 && u_coord < 399) && (v_coord > 1 && v_coord < 231))
        
        % Refer wikipedia of Bilinear interpolation
        % Q11 = (u_coord-1, v_coord-1), Q21 = (u_coord+1, v_coord-1), 
        % Q21 = (u_coord-1, v_coord+1), Q22 = (u_coord+1, v_coord+1)
       
        x1 = U_coord-1;
        x2 = U_coord+1;
        y1 = V_coord-1;
        y2 = V_coord+1;
        x = U_coord;
        y = V_coord;
            
        f_Q11 = color_warped(u_coord-1, v_coord-1,:);
        f_Q21 = color_warped(u_coord+1, v_coord-1,:);
        f_Q12 = color_warped(u_coord-1, v_coord+1,:);
        f_Q22 = color_warped(u_coord+1, v_coord+1,:);
        
        f_xy1 = ((x2-x)/(x2-x1))* f_Q11 + ((x-x1)/(x2-x1))* f_Q21;
        f_xy2 = ((x2-x)/(x2-x1))* f_Q12 + ((x-x1)/(x2-x1))* f_Q22;
        
        f_xy = ((y2-y)/(y2-y1))* f_xy1 + ((y-y1)/(y2-y1))* f_xy2;
        
      
        image_unwarped(i,j,:) = f_xy;
      end
    end
end
figure(3)
imshow(image_unwarped);
