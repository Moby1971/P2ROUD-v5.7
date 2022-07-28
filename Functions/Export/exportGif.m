function folderName = exportGif(app, tag)

% ------------------------------------
% ----- Export of images to gif ------
% ------------------------------------

% Proud data parameters object
obj = app.proudDataPars;

% Input images
gifImages = obj.images;

% Create folder if not exist, and delete folder content
folderName = [app.gifExportPath,[filesep,tag,'P-GIF']];
if (~exist(folderName, 'dir')); mkdir(folderName); end
delete([folderName,filesep,'*']);

% Phase orientation
if obj.PHASE_ORIENTATION == 1
    app.TextMessage('INFO: phase orientation = 1');
    gifImages = permute(rot90(permute(gifImages,[2 1 3 4 5 6]),1),[2 1 3 4 5 6]);
end

% Size of the data
dimx = size(gifImages,1);
dimy = size(gifImages,2);
dimz = size(gifImages,3);
dimd = size(gifImages,4);
NFA = size(gifImages,5);
NE = size(gifImages,6);

% Window, level, and scale between 0 and 255 grayvalues
window = app.WindowEditField.Value;
level = app.LevelEditField.Value;
window = window*255/max(gifImages(:));
level = level*255/max(gifImages(:));
gifImages = gifImages*255/max(gifImages(:));
gifImages = (255/window)*(gifImages - level + window/2);
gifImages(gifImages < 0) = 0;
gifImages(gifImages > 255) = 255;

% Aspect ratio
aspectRatio = (app.FOVViewField1.Value/app.FOVViewField2.Value)*(app.XEditField.Value/app.YEditField.Value);
if obj.PHASE_ORIENTATION
    aspectRatio = 1/aspectRatio;
end

% Oversample to increase image size, and aspect ratio
numRows = 2*dimx;
numCols = 2*round(dimy*aspectRatio);

% Export the gif images

for i=1:dimd % loop over all repetitions

    for z=1:dimz    % loop over all slices

        for j=1:NFA      % loop over all flip angles

            cine = false;
            if (obj.frame_loop_on == 1)
                cine = true;
            end

            if cine

                % File name
                fname = [folderName,filesep,'movie_d',num2str(i),'_s',num2str(z),'_fa',num2str(j),'.gif'];

                % Delay time
                delayTime = 1/NE;

                for k=1:NE      % loop over all echo times (cine loop)

                    % The image
                    im = rot90(uint8(round(imresize(squeeze(gifImages(:,:,z,i,j,k)),[numRows numCols]))));

                    % Write the gif file
                    if k==1
                        imwrite(im, fname,'DelayTime',delayTime,'LoopCount',inf);
                    else
                        imwrite(im, fname,'DelayTime',delayTime,'WriteMode','append','DelayTime',delayTime);
                    end

                end

            else

                % Delay time
                delayTime = 1/dimd;

                for k=1:NE      % loop over all echo times, multi-echo

                    % File name
                    fname = [folderName,filesep,'image_s',num2str(z),'_fa',num2str(j),'_te',num2str(k),'.gif'];

                    % The image
                    im = rot90(uint8(round(imresize(squeeze(gifImages(:,:,z,i,j,k)),[numRows numCols]))));

                    % Write the gif file
                    if i==1
                        imwrite(im, fname,'DelayTime',delayTime,'LoopCount',inf);
                    else
                        imwrite(im, fname,'DelayTime',delayTime,'WriteMode','append','DelayTime',delayTime);
                    end

                end

            end

        end

    end

end



% ------------------------------------
% ------- Flow images export ---------
% ------------------------------------

if obj.validFlow_flag

    % Input images
    gifImages = obj.flowImages;

    % Create folder if not exist, and delete folder content
    folderName = strcat(app.gifExportPath,filesep,'flowGIF');
    if (~exist(folderName, 'dir')); mkdir(folderName); end
    delete([folderName,filesep,'*']);

    % Phase orientation
    if obj.PHASE_ORIENTATION == 1
        gifImages = permute(rot90(permute(gifImages,[2 1 3 4 5 6 7]),1),[2 1 3 4 5 6 7]);
    end

    % Size of the data
    dimx = size(gifImages,1);
    dimy = size(gifImages,2);
    dimz = size(gifImages,3);
    NFA = size(gifImages,5);
    NE = size(gifImages,6);
    NF = size(gifImages,7);

    % Scale between 0 and 255 grayvalues
    gifImages(:,:,:,:,:,:,1) = round(127+gifImages(:,:,:,:,:,:,1)*128/abs(max(gifImages(:,:,:,:,:,:,1),[],'all')));
    gifImages(:,:,:,:,:,:,2) = round(127+gifImages(:,:,:,:,:,:,2)*128/abs(max(gifImages(:,:,:,:,:,:,2),[],'all')));
    gifImages(:,:,:,:,:,:,3) = round(127+gifImages(:,:,:,:,:,:,3)*128/abs(max(gifImages(:,:,:,:,:,:,3),[],'all')));
    gifImages(:,:,:,:,:,:,4) = round(127+gifImages(:,:,:,:,:,:,4)*128/abs(max(gifImages(:,:,:,:,:,:,4),[],'all')));
    gifImages(:,:,:,:,:,:,5) = round(    gifImages(:,:,:,:,:,:,5)*255/abs(max(gifImages(:,:,:,:,:,:,5),[],'all')));
    gifImages(gifImages < 0) = 0;
    gifImages(gifImages > 255) = 255;

    % Aspect ratio
    aspectRatio = (app.FOVViewField1.Value/app.FOVViewField2.Value)*(app.XEditField.Value/app.YEditField.Value);
    if obj.PHASE_ORIENTATION
        aspectRatio = 1/aspectRatio;
    end

    % Oversample to increase image size, and aspect ratio
    numRows = 2*dimx;
    numCols = 2*round(dimy*aspectRatio);

    % Export the gif images

    for i=1:NF % loop over flow directions

        for z=1:dimz    % loop over all slices

            for j=1:NFA      % loop over all flip angles

                % File name
                fname = [folderName,filesep,'movie_f',num2str(i),'_s',num2str(z),'_fa',num2str(j),'.gif'];

                % Delay time
                delayTime = 1/NE;

                for k=1:NE      % loop over all echo times (cine loop)

                    % The image
                    im = rot90(uint8(round(imresize(squeeze(gifImages(:,:,z,1,j,k,i)),[numRows numCols]))));

                    % Write the gif file
                    if k==1
                        imwrite(im, fname,'DelayTime',delayTime,'LoopCount',inf);
                    else
                        imwrite(im, fname,'DelayTime',delayTime,'WriteMode','append','DelayTime',delayTime);
                    end

                end

            end

        end

    end

end



end

