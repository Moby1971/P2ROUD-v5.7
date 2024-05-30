%% ---------------------------
%
%
%
%

clearvars;
close all;







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

