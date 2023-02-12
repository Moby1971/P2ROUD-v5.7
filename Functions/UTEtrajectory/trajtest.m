

traj = load('ktraj.txt');

%plot(traj)

gl = 128;

grad = ones(gl,1);

% 53 seems to overlap
pp = 55;

for i = 1:pp
    grad(i) = (i-1)/pp;
end

traj2 = zeros(gl,1);
for i = 1:gl
    traj2(i) = sum(grad(1:i));
end
traj2 = 0.5*traj2/max(traj2(:));

figure(1)
hold on

plot(traj)
plot(traj2)

hold off

filename = 'ktraj.txt';
fileID = fopen(filename,'w');

for i = 1:gl
   
    fprintf(fileID,num2str(traj2(i,1)));
    fprintf(fileID,'\n');
    
end

fclose(fileID);




