function folderName = exportDicomDCM(app, dcmdir)

% ---------------------------------
% Dicom export of images 
% Gustav Strijkers
% 31 March 2025
%
% ---------------------------------------------------------------


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
dimX = size(image,1);
dimY = size(image,2);
dimZ = size(image,3);
dimD = size(image,4);
dimFA = size(image,5);
dimE = size(image,6);
threeDdata = contains(lower(obj.dataType),lower(["3D","3Dp2roud","3Dute","ZTE"]));
seriesInstanceID = dicomuid;
     

% Reading in the DICOM header information
listing = dir(fullfile(dcmdir, '*.dcm'));
dcmFilename = [listing(1).folder,filesep,listing(1).name];
baseHeader = dicominfo(dcmFilename);
app.TextMessage(strcat('Reading DICOM info from',{' '},dcmFilename));

% Spatial z positions for 3D data
if threeDdata

    [~,spatial,~] = dicomreadVolume(listing(1).folder);
    spatialPos = unique(spatial.PatientPositions,'rows');

    dimzOrg = size(spatialPos,1);
    dimzRange = 1:(dimzOrg-1)/(dimZ-1):dimzOrg;
    
    v1 = interp1(spatialPos(:,1),dimzRange);
    v2 = interp1(spatialPos(:,2),dimzRange);
    v3 = interp1(spatialPos(:,3),dimzRange);

    spatialPositions = [v1',v2',v3'];

end

% Create new directory
ready = false;
cnt = 1;
while ~ready
    folderName = strcat(directory,filesep,"DICOM",filesep,num2str(baseHeader.SeriesNumber),'P',filesep,num2str(cnt),filesep);
    if ~exist(folderName, 'dir')
        mkdir(folderName);
        ready = true;
    end
    cnt = cnt + 1;
end

% export the dicom images
fileCounter = 0;
app.ExportProgressGauge.Value = 0;
totalNumberofImages = dimD*dimFA*dimE*dimZ;

for dynamic = 1:dimD                    % loop over all repetitions

    for flipAngle = 1:dimFA             % loop over all flip angles

        for echo = 1:dimE               % loop over all echo times

            for slice = 1:dimZ          % loop over all slices

                % Counter
                fileCounter = fileCounter + 1;

                % Read dicom header
             %   loc = dimZ+1-slice;
             %   dcmFilename = strcat(listing(loc).folder,filesep,listing(loc).name);
             %   baseHeader = dicominfo(dcmFilename);
             
                % File name
                fn = strcat('000000',num2str(fileCounter));
                fn = fn(size(fn,2)-5:size(fn,2));
                fname = strcat(folderName,filesep,fn,'.dcm');
                
                % Generate dicom header
                dcmHeader = generateDicomheaderDCM(baseHeader,fname,fileCounter,dynamic,flipAngle,echo,slice);

                % The image
                im = rot90(squeeze(cast(round(image(:,:,slice,dynamic,flipAngle,echo)),'uint16')));

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
    dimX = size(image,1);
    dimY = size(image,2);
    dimZ = size(image,3);
    dimD = 1;
    dimFA = size(image,5);
    dimE = size(image,6);
    dimFL = size(image,7);

    % For flow data the TE dimensions is the movie frame and TR the frame time
    app.pd.nr_frames = dimE;
    app.pd.timeperframe = obj.tr;

    % reading in the DICOM header information
    listing = dir(fullfile(dcmdir, '*.dcm'));
    dcmFilename = [listing(1).folder,filesep,listing(1).name];
    baseHeader = dicominfo(dcmFilename);
    app.TextMessage(strcat('Reading DICOM info from',{' '},dcmFilename));

    % create folder if not exist, and delete folder content
    folderName = stract(directory,filesep,"DICOM",filesep,num2str(baseHeader.SeriesNumber),'PF',num2str(dimD));
    if (~exist(folderName, 'dir'))
        mkdir(folderName);
    end
    delete(strcat(folderName,filesep,'*'));

    % export the dicom images
    fileCounter = 0;
    app.ExportProgressGauge.Value = 0;
    totalNumberofImages = dimFL*dimFA*dimE*dimZ;

    for dynamic=1:dimFL      % loop over all flow directions

        for flipAngle=1:dimFA      % loop over all flip angles

            for echo=1:dimE      % loop over all echo times

                for slice=1:dimZ        % loop over all slices

                    % Counter
                    fileCounter = fileCounter + 1;

                    % File name
                    fn = strcat('000000',num2str(fileCounter));
                    fn = fn(size(fn,2)-5:size(fn,2));
                    fname = strcat(folderName,filesep,'F',num2str(dynamic),fn,'.dcm');

                    % Dicom header
                    dcmHeader = generateDicomheaderDCMflow(baseHeader,fname,fileCounter,flipAngle,echo,slice);

                    % The image
                    im = rot90(squeeze(cast(round(image(:,:,slice,1,flipAngle,echo,dynamic)),'int16')));

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

    function dicomHeader = generateDicomheaderDCM(dcmHead,fname,~,dynamic,flipAngle,echo,slice)

        % GENERATES DICOM HEADER FOR EXPORT
        %
        % dcmHead = dicom info from scanner generated dicom

        frametime = 1000*obj.acqdur/dimD;    % time between frames in ms

        if obj.PHASE_ORIENTATION == 1
            pixely = app.FOVViewField1.Value/dimY;
            pixelx = app.FOVViewField2.Value/dimX;
        else
            pixely = app.FOVViewField2.Value/dimY;
            pixelx = app.FOVViewField1.Value/dimX;
        end

        sliceSep = obj.SLICE_SEPARATION + obj.SLICE_THICKNESS;

        dcmHead.RepetitionTime = obj.TR;
        dcmHead.Filename = fname;
        dcmHead.FileModDate = obj.date;
        dcmHead.FileSize = dimY*dimX*2;
        dcmHead.Width = dimY;
        dcmHead.Height = dimX;
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
        dcmHead.TriggerTime = (dynamic-1)*frametime;    % frame time (ms)
        dcmHead.AcquisitionMatrix = uint16([dimX 0 0 dimY])';
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.AcquisitionDuration = obj.acqdur;
        % dcmHead.InstanceNumber = fileCounter;          
        dcmHead.InstanceNumber = (dynamic-1)*dimZ*dimE + (echo-1)*dimZ + slice;     % instance number
            
        dcmHead.TemporalPositionIdentifier = dynamic;
        dcmHead.NumberOfTemporalPositions = dimD;
        dcmHead.ImagesInAcquisition = dimD*dimZ*dimE;
        % dcmHead.TemporalPositionIndex = dynamic;
        dcmHead.TemporalPositionIndex = (dynamic-1)*dimE + echo;
           
        dcmHead.Rows = dimY;
        dcmHead.Columns = dimX;
        dcmHead.PixelSpacing = [pixely pixelx]';
        dcmHead.PixelAspectRatio = [1 pixely/pixelx]';
        dcmHead.BitsAllocated = 16;
        dcmHead.BitsStored = 15;
        dcmHead.HighBit = 14;
        dcmHead.PixelRepresentation = 0;
        dcmHead.PixelPaddingValue = 0;
        dcmHead.RescaleIntercept = 0;
        dcmHead.RescaleSlope = 1;
        dcmHead.NumberOfSlices = dimZ;

        dcmHead.SliceThickness = obj.SLICE_THICKNESS;
        if threeDdata
            dcmHead.SpacingBetweenSlices = obj.SLICE_THICKNESS;
        else
            dcmHead.SpacingBetweenSlices = sliceSep/obj.SLICE_INTERLEAVE;
        end

        dcmHead.EchoTrainLength = obj.NO_ECHOES;
        dcmHead.EchoTime = obj.TE*echo;                         % ECHO TIME
        dcmHead.FlipAngle = obj.flipAngleArray(flipAngle);           % FLIP ANGLES

        if threeDdata

            dcmHead.ImagePositionPatient(1) = spatialPositions(slice,1);
            dcmHead.ImagePositionPatient(2) = spatialPositions(slice,2);
            dcmHead.ImagePositionPatient(3) = spatialPositions(slice,3);

        else

           startslice = 0;
           if isfield(dcmHead, 'SliceLocation')
               startslice = dcmHead.SliceLocation;
           end
           dcmHead.SliceLocation = startslice+(slice-1)*(sliceSep/obj.SLICE_INTERLEAVE);

        end

        dcmHead.SeriesInstanceUID = seriesInstanceID;
      
        dicomHeader = dcmHead;

    end % Generate dicom header


    function dicomHeader = generateDicomheaderDCMflow(dcmHead,fname,fileCounter,j,k,z)

        % GENERATES DICOM HEADER FOR EXPORT OF FLOW IMAGES
        %
        % dcmHead = dicom info from scanner generated dicom

        frametime = 1000*obj.acqdur/dimD;    % time between frames in ms

        if obj.PHASE_ORIENTATION == 1
            pixely = app.FOVViewField1.Value/dimY;
            pixelx = app.FOVViewField2.Value/dimX;
        else
            pixely = app.FOVViewField2.Value/dimY;
            pixelx = app.FOVViewField1.Value/dimX;
        end

        threeDdata = contains(obj.dataType,["3D","3Dp2roud","3Dute"]);
        sliceSep = obj.SLICE_SEPARATION + obj.SLICE_THICKNESS;

        dcmHead.RepetitionTime = obj.TR;
        dcmHead.Filename = fname;
        dcmHead.FileModDate = obj.date;
        dcmHead.FileSize = dimY*dimX*2;
        dcmHead.Width = dimY;
        dcmHead.Height = dimX;
        dcmHead.BitDepth = 15;
        dcmHead.InstitutionName = 'Amsterdam UMC - MRI';
        dcmHead.ReferringPhysicianName.FamilyName = 'Amsterdam UMC - MRI';
        dcmHead.StudyDescription = 'CINE';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC - MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC - MRI';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC - MRI';
        dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC - MRI';
        dcmHead.OperatorName.FamilyName = 'Amsterdam UMC - MRI';
        dcmHead.ManufacturerModelName = 'MR Solutions';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.NumberOfAverages = obj.NO_AVERAGES;
        dcmHead.InversionTime = 0;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.TriggerTime = (k-1)*frametime;                      % frame time (ms)
        dcmHead.AcquisitionMatrix = uint16([dimX 0 0 dimY])';
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.AcquisitionDuration = obj.acqdur;
        dcmHead.InstanceNumber = fileCounter;                       % instance number
        dcmHead.TemporalPositionIdentifier = k;                     % frame number
        dcmHead.NumberOfTemporalPositions = dimD;
        dcmHead.ImagesInAcquisition = dimD*dimZ;
        dcmHead.TemporalPositionIndex = k;
        dcmHead.Rows = dimY;
        dcmHead.Columns = dimX;
        dcmHead.PixelSpacing = [pixely pixelx]';
        dcmHead.PixelAspectRatio = [1 pixely/pixelx]';
        dcmHead.BitsAllocated = 16;
        dcmHead.BitsStored = 15;
        dcmHead.HighBit = 14;
        dcmHead.PixelRepresentation = 0;
        dcmHead.PixelPaddingValue = 0;
        dcmHead.RescaleIntercept = 0;
        dcmHead.RescaleSlope = 1;
        dcmHead.NumberOfSlices = dimZ;
       
        dcmHead.EchoTime = obj.TE*k;                                                    % ECHO TIME
        dcmHead.EchoTrainLength = obj.NO_ECHOES;
        dcmHead.FlipAngle = obj.flipAngleArray(j);                                      % FLIP ANGLES

        dcmHead.SliceThickness = obj.SLICE_THICKNESS;
        if threeDdata
            dcmHead.SpacingBetweenSlices = obj.SLICE_THICKNESS;
        else
            dcmHead.SpacingBetweenSlices = sliceSep/obj.SLICE_INTERLEAVE;
        end
    
        if threeDdata
            dcmHead.ImagePositionPatient(2) = dcmHead.ImagePositionPatient(2) + (z-1)*obj.SLICE_THICKNESS;
            dcmHead.SliceLocation = dcmHead.ImagePositionPatient(2);
        else
            startslice = 0;
            if isfield(dcmHead, 'SliceLocation')
                startslice = dcmHead.SliceLocation;
            end
            dcmHead.SliceLocation = startslice+(z-1)*(sliceSep/obj.SLICE_INTERLEAVE);
        end

        dcmHead.SeriesInstanceUID = seriesInstanceID;
        dcmHead.SequenceVariant = 'NONE';
        dcmHead.ScanOptions = 'CG';
        dcmHead.MRAcquisitionType = '2D';

        dicomHeader = dcmHead;

    end % Generate dicom header flow

% ---- END OF FUNCTIONS ------


end % exportDicomDCM