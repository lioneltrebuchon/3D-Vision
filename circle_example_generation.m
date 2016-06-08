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
imwrite(cir,'/Users/markopanjek/Documents/DOKUMENTI/Sola/ETH/3D_Vision/semester_project/circle_original.jpeg');



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
    i
end
figure(2)
imshow(image_final);
imwrite(image_final,'/Users/markopanjek/Documents/DOKUMENTI/Sola/ETH/3D_Vision/semester_project/circle_warped.jpeg');
            

%% circle camera view
clear image_final;
theta = 0;
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
    i
end
figure(2)
imshow(image_final);
%imwrite(image_final,'/Users/markopanjek/Documents/DOKUMENTI/Sola/ETH/3D_Vision/semester_project/circle_warped.jpeg');
            
%% UNWARPING USING HOMOGRAPHY FROM C++ (Torstens meeting, 8.6.2016)
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
        %if((u_coord > 0 && u_coord < 400) && (v_coord > 0 && v_coord < 400))
        image_unwarped(i,j) = image_final(u_coord,v_coord);
        %end
        
    end
end
figure(4)
imshow(image_unwarped);