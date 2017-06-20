%Simple demos for the FastICP code
%Use the classic point cloud Bunny Rabbit created by Stanford 3D Scanning
%Repository to demo

%Xu CHEN <cx495086@outlook.com>
%2017.06

system('FastICP.exe');

load 'model.txt'
load 'data.txt';
load 'RotateT.txt';

modelmatched=model*RotateT(1:3,1:3)'+ones(size(model,1),1)*RotateT(1:3,4)';

%% show
figure
showPointCloud(model,[1 0.5 0],'MarkerSize' ,2);
hold on;
showPointCloud(data,[0 0.5 1],'MarkerSize' ,2);
showPointCloud(modelmatched,[1 0 0],'MarkerSize' ,2);

xlabel('X');ylabel('Y');zlabel('Z');
axis equal;
rotate3d on;
grid off;
