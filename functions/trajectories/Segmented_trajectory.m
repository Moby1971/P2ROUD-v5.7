%% ---------------------------
%
%
%
%

clearvars


traj = zeros(8*12,1);

for j = 1:64
    l(j) = j-33;
end
l(l==-2) = [];
l(l==-1) = [];
l(l==0) = [];
l(l==1) = [];


for i = 1:8
    traj(1+(i-1)*12) = l(randi(60));
    traj(2+(i-1)*12) = l(randi(60));
    traj(3+(i-1)*12) = l(randi(60));
    traj(4+(i-1)*12) = -2;
    traj(5+(i-1)*12) = -1;
    traj(6+(i-1)*12) = 0;
    traj(7+(i-1)*12) = 1;
    traj(8+(i-1)*12) = l(randi(60));
    traj(9+(i-1)*12) =  l(randi(60));
    traj(10+(i-1)*12) = l(randi(60));
    traj(11+(i-1)*12) = l(randi(60));
    traj(12+(i-1)*12) = l(randi(60));
end


for i = 1:8
    
    
    
   trajs((i-1)*12+1 : (i-1)*12 + 12) = sort(traj((i-1)*12+1 : (i-1)*12 + 12),'ascend');
   
   if rem(i, 2)== 0 
       trajs((i-1)*12+1 : (i-1)*12 + 12) = sort(traj((i-1)*12+1 : (i-1)*12 + 12),'descend');
   end
   
   
end

figure(1);
plot(trajs);


kspace = zeros(11,64,64);

cnt = 1;

for i = 1:8
 for j = 1:12
    kspace(i,trajs(cnt)+33,:) = 1;
    cnt = cnt + 1;
 end
 
end

kspace(9,:,:) = kspace(1,:,:) + kspace(2,:,:) + kspace(3,:,:) + kspace(4,:,:) + kspace(5,:,:) + kspace(6,:,:) + kspace(7,:,:) + kspace(8,:,:);
kspace(10,:,:) = kspace(1,:,:) + kspace(2,:,:) + kspace(3,:,:) + kspace(4,:,:);
kspace(11,:,:) = kspace(5,:,:) + kspace(6,:,:) + kspace(7,:,:) + kspace(8,:,:);


ptitle = ['1','2','3','4','5','6','7','8',"1-8",'1-4','5-8'];


figure(2);

for i = 1:11
    subplot(2,6,i),imshow(squeeze(kspace(i,:,:)),[0 1]);
    title(ptitle(i));
end


filename = 'per_dimy96_segm8_lines12.txt';
fileID = fopen(filename,'w');

for i = 1:length(trajs)
   
    fprintf(fileID,num2str(trajs(1,i)));
    fprintf(fileID,'\n');
    
end

fclose(fileID);

