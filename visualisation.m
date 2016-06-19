file=fopen('ETHday.csv');
all = textscan(file,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Headerlines',1,'delimiter',',');
fclose(file);
all = cell2mat(all);
points = all(:,1:3);
normals = all(:,4:6);
camera1s = all(:,7:9);
camera2s = all(:,10:12);
RWto1s = all(:,13:21);
R2toWs = all(:,22:30);
numPoints = size(all,1);
for i = 1:numPoints
    RWto1 = reshape(RWto1s(i,:),[3 3])';
    R2toW = reshape(R2toWs(i,:),[3 3])';
    vector_visualizations(R2toW, RWto1,camera1s(i,:),camera2s(i,:),normals(i,:),points(i,:));
end
% skip = 50;
% figure
% quiver3(points(1:skip:end,1),points(1:skip:end,2),points(1:skip:end,3),normals(1:skip:end,1),normals(1:skip:end,2),normals(1:skip:end,3));
% view(2);
% figure
% scatter3(points(1:skip:end,1),points(1:skip:end,2),points(1:skip:end,3));
% view(2);