function folderName = exportDicomDCM(app, dcmdir)

% ---------------------------------
% ---- Dicom export of images -----
% ---------------------------------


% Proud data parameters object
obj = app.pd;

directory = app.dicomExportPath;
image = obj.images;

% Phase orientation
if obj.PHASE_ORIENTATION == 1
    app.TextMessage('INFO: phase orientation = 1');
    image = permute(rot90(permute(image,[2 1 3 4 5 6]),1),[2 1 3 4 5 6]);
end

% size of the data
dimx = size(image,1);
dimy = size(image,2);
dimz = size(image,3);
dimd = size(image,4);
NFA = size(image,5);
NE = size(image,6);

% reading in the DICOM header information
listing = dir(fullfile(dcmdir, '*.dcm'));
dcmFilename = [listing(1).folder,filesep,listing(1).name];
baseHeader = dicominfo(dcmFilename);
app.TextMessage(strcat('Reading DICOM info from',{' '},dcmFilename));

% Create new directory
ready = false;
cnt = 1;
while ~ready
    folderName = strcat(directory,filesep,"DICOM",filesep,num2str(baseHeader.SeriesNumber),'P',num2str(dimd),filesep,num2str(cnt),filesep);
    if ~exist(folderName, 'dir')
        mkdir(folderName);
        ready = true;
    end
    cnt = cnt + 1;
end

% export the dicom images
fileCounter = 0;
app.ExportProgressGauge.Value = 0;
totalNumberofImages = dimd*NFA*NE*dimz;

for i=1:dimd % loop over all repetitions

    for j=1:NFA      % loop over all flip angles

        for k=1:NE      % loop over all echo times

            for z=1:dimz        % loop over all slices

                % Counter
                fileCounter = fileCounter + 1;

                % File name
                fn = strcat('000000',num2str(fileCounter));
                fn = fn(size(fn,2)-5:size(fn,2));
                fname = strcat(folderName,filesep,fn,'.dcm');

                % Dicom header
                dcmHeader = generateDicomheaderDCM(baseHeader,fname,fileCounter,i,j,k,z);

                % The image
                im = rot90(squeeze(cast(round(image(:,:,z,i,j,k)),'uint16')));

                % Write the dicom file
                dicomwrite(im, fname, dcmHeader);

                % Update progress bar
                app.ExportProgressGauge.Value = round(100*fileCounter/totalNumberofImages);
                drawnow;

            end

        end

    end

end


% --------------------------------------
% ---- Dicom export of flow images -----
% --------------------------------------

if obj.validFlow_flag

    app.TextMessage('Exporting flow images to dicom ...');

    % Proud data parameters object
    obj = app.pd;

    directory = app.dicomExportPath;
    image = obj.flowImages*100;  % scale by factor 100

    % Phase orientation
    if obj.PHASE_ORIENTATION == 1
        image = permute(rot90(permute(image,[2 1 3 4 5 6 7]),1),[2 1 3 4 5 6 7]);
    end

    % size of the data
    dimx = size(image,1);
    dimy = size(image,2);
    dimz = size(image,3);
    dimd = 1;
    NFA = size(image,5);
    NE = size(image,6);
    NF = size(image,7);

    % For flow data the TE dimensions is the movie frame and TR the frame time
    app.pd.nr_frames = NE;
    app.pd.timeperframe = app.pd.tr;

    % reading in the DICOM header information
    listing = dir(fullfile(dcmdir, '*.dcm'));
    dcmFilename = [listing(1).folder,filesep,listing(1).name];
    baseHeader = dicominfo(dcmFilename);
    app.TextMessage(strcat('Reading DICOM info from',{' '},dcmFilename));

    % create folder if not exist, and delete folder content
    folderName = stract(directory,filesep,"DICOM",filesep,num2str(baseHeader.SeriesNumber),'PF',num2str(dimd));
    if (~exist(folderName, 'dir'))
        mkdir(folderName);
    end
    delete(strcat(folderName,filesep,'*'));

    % export the dicom images
    fileCounter = 0;
    app.ExportProgressGauge.Value = 0;
    totalNumberofImages = NF*NFA*NE*dimz;

    for i=1:NF      % loop over all flow directions

        for j=1:NFA      % loop over all flip angles

            for k=1:NE      % loop over all echo times

                for z=1:dimz        % loop over all slices

                    % Counter
                    fileCounter = fileCounter + 1;

                    % File name
                    fn = strcat('000000',num2str(fileCounter));
                    fn = fn(size(fn,2)-5:size(fn,2));
                    fname = strcat(folderName,filesep,'F',num2str(i),fn,'.dcm');

                    % Dicom header
                    dcmHeader = generateDicomheaderDCMflow(baseHeader,fname,fileCounter,j,k,z);

                    % The image
                    im = rot90(squeeze(cast(round(image(:,:,z,1,j,k,i)),'int16')));

                    % Write the dicom file
                    dicomwrite(im, fname, dcmHeader);

                    % Update progress bar
                    app.ExportProgressGauge.Value = round(100*fileCounter/totalNumberofImages);
                    drawnow;

                end

            end

        end

    end

end



% ----------------------
% --- FUNCTIONS --------
% ----------------------

    function dicomHeader = generateDicomheaderDCM(dcmHead,fname,fileCounter,i,j,k,z)

        % GENERATES DICOM HEADER FOR EXPORT
        %
        % dcmHead = dicom info from scanner generated dicom

        frametime = 1000*obj.acqdur/dimd;    % time between frames in ms

        if obj.PHASE_ORIENTATION == 1
            pixely = app.FOVViewField1.Value/dimy;
            pixelx = app.FOVViewField2.Value/dimx;
        else
            pixely = app.FOVViewField2.Value/dimy;
            pixelx = app.FOVViewField1.Value/dimx;
        end

        dcmHead.RepetitionTime = obj.TR;
        dcmHead.Filename = fname;
        dcmHead.FileModDate = obj.date;
        dcmHead.FileSize = dimy*dimx*2;
        dcmHead.Width = dimy;
        dcmHead.Height = dimx;
        dcmHead.BitDepth = 15;
        dcmHead.InstitutionName = 'Amsterdam UMC';
        dcmHead.ReferringPhysicianName.FamilyName = 'AMC preclinical MRI';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.OperatorName.FamilyName = 'manager';
        dcmHead.ManufacturerModelName = 'MRS7024';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.NumberOfAverages = obj.NO_AVERAGES;
        dcmHead.InversionTime = 0;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.MagneticFieldStrength = 7;
        dcmHead.TriggerTime = (i-1)*frametime;    % frame time (ms)
        dcmHead.AcquisitionMatrix = uint16([dimx 0 0 dimy])';
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.AcquisitionDuration = obj.acqdur;
        dcmHead.InstanceNumber = fileCounter;          % instance number
        dcmHead.TemporalPositionIdentifier = i;     % frame number
        dcmHead.NumberOfTemporalPositions = dimd;
        dcmHead.ImagesInAcquisition = dimd*dimz;
        dcmHead.TemporalPositionIndex = i;
        dcmHead.Rows = dimy;
        dcmHead.Columns = dimx;
        dcmHead.PixelSpacing = [pixely pixelx]';
        dcmHead.PixelAspectRatio = [1 pixely/pixelx]';
        dcmHead.BitsAllocated = 16;
        dcmHead.BitsStored = 15;
        dcmHead.HighBit = 14;
        dcmHead.PixelRepresentation = 0;
        dcmHead.PixelPaddingValue = 0;
        dcmHead.RescaleIntercept = 0;
        dcmHead.RescaleSlope = 1;
        dcmHead.NumberOfSlices = dimz;

        dcmHead.SliceThickness = obj.SLICE_THICKNESS;
        dcmHead.EchoTime = obj.TE*k;                 % ECHO TIME
        dcmHead.SpacingBetweenSlices = obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE;
        dcmHead.EchoTrainLength = obj.NO_ECHOES;
        dcmHead.FlipAngle = obj.flipAngleArray(j);           % FLIP ANGLES

        if isfield(dcmHead, 'SliceLocation')
            startslice = dcmHead.SliceLocation;
            dcmHead.SliceLocation = startslice+(z-1)*(obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE);
        end

        dicomHeader = dcmHead;

    end % Generate dicom header


    function dicomHeader = generateDicomheaderDCMflow(dcmHead,fname,fileCounter,j,k,z)

        % GENERATES DICOM HEADER FOR EXPORT OF FLOW IMAGES
        %
        % dcmHead = dicom info from scanner generated dicom

        frametime = 1000*obj.acqdur/dimd;    % time between frames in ms

        if obj.PHASE_ORIENTATION == 1
            pixely = app.FOVViewField1.Value/dimy;
            pixelx = app.FOVViewField2.Value/dimx;
        else
            pixely = app.FOVViewField2.Value/dimy;
            pixelx = app.FOVViewField1.Value/dimx;
        end

        dcmHead.RepetitionTime = obj.TR;
        dcmHead.Filename = fname;
        dcmHead.FileModDate = obj.date;
        dcmHead.FileSize = dimy*dimx*2;
        dcmHead.Width = dimy;
        dcmHead.Height = dimx;
        dcmHead.BitDepth = 15;
        dcmHead.InstitutionName = 'Amsterdam UMC';
        dcmHead.ReferringPhysicianName.FamilyName = 'AMC preclinical MRI';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.OperatorName.FamilyName = 'manager';
        dcmHead.ManufacturerModelName = 'MRS7024';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.NumberOfAverages = obj.NO_AVERAGES;
        dcmHead.InversionTime = 0;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.MagneticFieldStrength = 7;
        dcmHead.TriggerTime = (k-1)*frametime;    % frame time (ms)
        dcmHead.AcquisitionMatrix = uint16([dimx 0 0 dimy])';
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.AcquisitionDuration = obj.acqdur;
        dcmHead.InstanceNumber = fileCounter;          % instance number
        dcmHead.TemporalPositionIdentifier = k;     % frame number
        dcmHead.NumberOfTemporalPositions = dimd;
        dcmHead.ImagesInAcquisition = dimd*dimz;
        dcmHead.TemporalPositionIndex = k;
        dcmHead.Rows = dimy;
        dcmHead.Columns = dimx;
        dcmHead.PixelSpacing = [pixely pixelx]';
        dcmHead.PixelAspectRatio = [1 pixely/pixelx]';
        dcmHead.BitsAllocated = 16;
        dcmHead.BitsStored = 15;
        dcmHead.HighBit = 14;
        dcmHead.PixelRepresentation = 0;
        dcmHead.PixelPaddingValue = 0;
        dcmHead.RescaleIntercept = 0;
        dcmHead.RescaleSlope = 1;
        dcmHead.NumberOfSlices = dimz;

        dcmHead.SliceThickness = obj.SLICE_THICKNESS;
        dcmHead.EchoTime = obj.TE*k;                 % ECHO TIME
        dcmHead.SpacingBetweenSlices = obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE;
        dcmHead.EchoTrainLength = obj.NO_ECHOES;
        dcmHead.FlipAngle = obj.flipAngleArray(j);           % FLIP ANGLES

        if isfield(dcmHead, 'SliceLocation')
            startslice = dcmHead.SliceLocation;
            dcmHead.SliceLocation = startslice+(z-1)*(obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE);
        end

        dicomHeader = dcmHead;

    end % Generate dicom header flow

% ---- END OF FUNCTIONS ------


end % exportDicomDCM