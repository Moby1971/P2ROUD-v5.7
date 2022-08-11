function folderName = exportDicomMRD(app, tag)



% -----------------------------------
% ----- DICOM export of images ------
% -----------------------------------


% Proud data parameters object
obj = app.proudDataPars;

directory = app.dicomExportPath;
image = obj.images;

% Create folder if not exist, and clear
folderName = [directory,[filesep,tag,'P']];
if (~exist(folderName, 'dir')); mkdir(folderName); end
delete([folderName,filesep,'*']);

% Phase orientation
if obj.PHASE_ORIENTATION == 1
    app.TextMessage('INFO: phase orientation = 1');
    image = permute(rot90(permute(image,[2 1 3 4 5 6]),1),[2 1 3 4 5 6]);
end

[dimx,dimy,dimz,NR,NFA,NE] = size(image);

% Export the dicom images
dcmid = dicomuid;   % unique identifier
dcmid = dcmid(1:50);

fileCounter = 0;
app.ExportProgressGauge.Value = 0;
totalNumberofImages = NR*NFA*NE*dimz;

for i = 1:NR      % loop over all repetitions

    for j = 1:NFA     % loop over all flip angles

        for k = 1:NE      % loop over all echo times

            for z = 1:dimz        % loop over all slices

                % Counter
                fileCounter = fileCounter + 1;

                % File name
                fn = ['00000',num2str(fileCounter)];
                fn = fn(size(fn,2)-5:size(fn,2));
                fname = [folderName,filesep,fn,'.dcm'];

                % Dicom header
                dcmHeader = generateDicomheaderMRD(fname,fileCounter,i,j,k,z);

                % The image
                image1 = rot90(squeeze(cast(round(image(:,:,z,i,j,k)),'uint16')));

                % Write the dicom file
                dicomwrite(image1, fname, dcmHeader);

                % Update progress bar
                app.ExportProgressGauge.Value = round(100*fileCounter/totalNumberofImages);
                drawnow;

            end

        end

    end

end




% ----------------------------------------
% ----- DICOM export of flow images ------
% ----------------------------------------


if obj.validFlow_flag

    app.TextMessage('Exporting flow images to dicom ...');

    % Proud data parameters object
    obj = app.proudDataPars;

    directory = app.dicomExportPath;
    image = obj.flowImages*100; %% scale by factor 100

    % Create folder if not exist, and clear
    folderName = [directory,[filesep,tag,'PF']];
    if (~exist(folderName, 'dir')); mkdir(folderName); end
    delete([folderName,filesep,'*']);

    % Phase orientation
    if obj.PHASE_ORIENTATION == 1
        image = permute(rot90(permute(image,[2 1 3 4 5 6 7]),1),[2 1 3 4 5 6 7]);
    end

    [dimx,dimy,dimz,~,NFA,NE,NF] = size(image);
    NR = 1;

    % Export the dicom images
    dcmid = dicomuid;   % unique identifier
    dcmid = dcmid(1:50);

    fileCounter = 0;
    app.ExportProgressGauge.Value = 0;
    totalNumberofImages = NF*NFA*NE*dimz;

    for i = 1:NF      % loop over all flow directions

        for j = 1:NFA     % loop over all flip angles

            for k = 1:NE      % loop over all echo times

                for z = 1:dimz        % loop over all slices

                    % Counter
                    fileCounter = fileCounter + 1;

                    % File name
                    fn = ['00000',num2str(fileCounter)];
                    fn = fn(size(fn,2)-5:size(fn,2));
                    fname = [folderName,filesep,'F',num2str(i),fn,'.dcm'];

                    % Dicom header
                    dcmHeader = generateDicomheaderMRDflow(fname,fileCounter,j,k,z);

                    % The image
                    image1 = rot90(squeeze(cast(round(image(:,:,z,1,j,k,i)),'int16')));

                    % Write the dicom file
                    dicomwrite(image1, fname, dcmHeader);

                    % Update progress bar
                    app.ExportProgressGauge.Value = round(100*fileCounter/totalNumberofImages);
                    drawnow;

                end

            end

        end

    end

end


app.ExportProgressGauge.Value = 100;



% ---------------------
% ---- functions ------
% ---------------------

    function dicomHeader = generateDicomheaderMRD(fname,fileCounter,i,j,k,z)

        % GENERATES DICOM HEADER FOR EXPORT

        acq_dur = obj.nr_frames * obj.timeperframe;   % acquisition time in seconds

        if obj.PHASE_ORIENTATION == 1
            pixelx = obj.aspectratio*obj.FOV/dimx;
            pixely = obj.FOV/dimy;
        else
            pixelx = obj.FOV/dimx;
            pixely = obj.aspectratio*obj.FOV/dimy;
        end

        dt = datetime(obj.date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
        year = num2str(dt.Year);
        month = ['0',num2str(dt.Month)]; month = month(end-1:end);
        day = ['0',num2str(dt.Day)]; day = day(end-1:end);
        date = [year,month,day];

        hour = ['0',num2str(dt.Hour)]; hour = hour(end-1:end);
        minute = ['0',num2str(dt.Minute)]; minute = minute(end-1:end);
        seconds = ['0',num2str(dt.Second)]; seconds = seconds(end-1:end);
        time = [hour,minute,seconds];

        dcmHead.Filename = fname;
        dcmHead.FileModDate = obj.date;
        dcmHead.FileSize = dimy*dimz*2;
        dcmHead.Format = 'DICOM';
        dcmHead.FormatVersion = 3;
        dcmHead.Width = dimx;
        dcmHead.Height = dimy;
        dcmHead.BitDepth = 15;
        dcmHead.ColorType = 'grayscale';
        dcmHead.FileMetaInformationGroupLength = 178;
        dcmHead.FileMetaInformationVersion = uint8([0, 1])';
        dcmHead.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.TransferSyntaxUID = '1.2.840.10008.1.2.1';
        dcmHead.ImplementationClassUID = '1.2.826.0.9717382.3.0.3.6.0';
        dcmHead.ImplementationVersionName = 'OFFIS_DCMTK_360';
        dcmHead.SpecificCharacterSet = 'ISO_IR 100';
        dcmHead.ImageType = 'DERIVED\4D-MRI\';
        dcmHead.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.StudyDate = date;
        dcmHead.SeriesDate = date;
        dcmHead.AcquisitionDate = date;
        dcmHead.StudyTime = time;
        dcmHead.SeriesTime = time;
        dcmHead.AcquisitionTime = 1000*(i-1)*obj.timeperframe;
        dcmHead.ContentTime = time;
        dcmHead.Modality = 'MR';
        dcmHead.Manufacturer = 'MR Solutions Ltd';
        dcmHead.InstitutionName = 'Amsterdam UMC';
        dcmHead.InstitutionAddress = 'Amsterdam, Netherlands';
        dcmHead.ReferringPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.ReferringPhysicianName.GivenName = '';
        dcmHead.ReferringPhysicianName.MiddleName = '';
        dcmHead.ReferringPhysicianName.NamePrefix = '';
        dcmHead.ReferringPhysicianName.NameSuffix = '';
        dcmHead.StationName = 'MRI Scanner';
        dcmHead.StudyDescription = 'XD-data';
        dcmHead.SeriesDescription = '';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.GivenName = '';
        dcmHead.PhysicianOfRecord.MiddleName = '';
        dcmHead.PhysicianOfRecord.NamePrefix = '';
        dcmHead.PhysicianOfRecord.NameSuffix = '';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PerformingPhysicianName.GivenName = '';
        dcmHead.PerformingPhysicianName.MiddleName = '';
        dcmHead.PerformingPhysicianName.NamePrefix = '';
        dcmHead.PerformingPhysicianName.NameSuffix = '';
        dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianReadingStudy.GivenName = '';
        dcmHead.PhysicianReadingStudy.MiddleName = '';
        dcmHead.PhysicianReadingStudy.NamePrefix = '';
        dcmHead.PhysicianReadingStudy.NameSuffix = '';
        dcmHead.OperatorName.FamilyName = 'manager';
        dcmHead.AdmittingDiagnosesDescription = '';
        dcmHead.ManufacturerModelName = 'MRS7024';
        dcmHead.ReferencedSOPClassUID = '';
        dcmHead.ReferencedSOPInstanceUID = '';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.DerivationDescription = '';
        dcmHead.FrameType = '';
        dcmHead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PatientID = '01';
        dcmHead.PatientBirthDate = date;
        dcmHead.PatientBirthTime = '';
        dcmHead.PatientSex = 'F';
        dcmHead.OtherPatientID = '';
        dcmHead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.OtherPatientName.GivenName = '';
        dcmHead.OtherPatientName.MiddleName = '';
        dcmHead.OtherPatientName.NamePrefix = '';
        dcmHead.OtherPatientName.NameSuffix = '';
        dcmHead.PatientAge = '1';
        dcmHead.PatientSize = [];
        dcmHead.PatientWeight = 0.0300;
        dcmHead.Occupation = '';
        dcmHead.AdditionalPatientHistory = '';
        dcmHead.PatientComments = '';
        dcmHead.BodyPartExamined = '';
        dcmHead.SequenceName = obj.PPL;
        dcmHead.SliceThickness = obj.SLICE_THICKNESS;
        dcmHead.KVP = 0;
        dcmHead.RepetitionTime = obj.TR;
        dcmHead.EchoTime = obj.TE*k;                 % ECHO TIME
        dcmHead.InversionTime = 0;
        dcmHead.NumberOfAverages = obj.NO_AVERAGES;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.MagneticFieldStrength = obj.field_strength;
        dcmHead.SpacingBetweenSlices = obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE;
        dcmHead.EchoTrainLength = obj.NO_ECHOES;
        dcmHead.DeviceSerialNumber = '0034';
        dcmHead.PlateID = '';
        dcmHead.SoftwareVersion = '1.0.0.0';
        dcmHead.ProtocolName = '';
        dcmHead.SpatialResolution = [];
        dcmHead.TriggerTime = 1000*(i-1)*obj.timeperframe;  % trigger time in ms
        dcmHead.DistanceSourceToDetector = [];
        dcmHead.DistanceSourceToPatient = [];
        dcmHead.FieldofViewDimensions = [obj.FOV obj.aspectratio*obj.FOV obj.SLICE_THICKNESS];
        dcmHead.ExposureTime = [];
        dcmHead.XrayTubeCurrent = [];
        dcmHead.Exposure = [];
        dcmHead.ExposureInuAs = [];
        dcmHead.FilterType = '';
        dcmHead.GeneratorPower = [];
        dcmHead.CollimatorGridName = '';
        dcmHead.FocalSpot = [];
        dcmHead.DateOfLastCalibration = '';
        dcmHead.TimeOfLastCalibration = '';
        dcmHead.PlateType = '';
        dcmHead.PhosphorType = '';
        dcmHead.AcquisitionMatrix = uint16([dimx 0 0 dimy])';
        dcmHead.FlipAngle = obj.flipAngleArray(j);           % FLIP ANGLES
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.CassetteOrientation = 'PORTRAIT';
        dcmHead.CassetteSize = '25CMX25CM';
        dcmHead.ExposuresOnPlate = 0;
        dcmHead.RelativeXrayExposure = [];
        dcmHead.AcquisitionComments = '';
        dcmHead.PatientPosition = 'HFS';
        dcmHead.Sensitivity = [];
        dcmHead.FieldOfViewOrigin = [];
        dcmHead.FieldOfViewRotation = [];
        dcmHead.AcquisitionDuration = acq_dur;
        dcmHead.StudyInstanceUID = dcmid(1:18);
        dcmHead.SeriesInstanceUID = [dcmid(1:18),'.',num2str(obj.filename)];
        dcmHead.StudyID = '01';
        dcmHead.SeriesNumber = obj.filename;
        dcmHead.AcquisitionNumber = 1;
        dcmHead.InstanceNumber = fileCounter;
        dcmHead.ImagePositionPatient = [-obj.FOV/2 -(obj.aspectratio*obj.FOV/2) (z-round(obj.NO_SLICES/2))*(obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE)]';
        dcmHead.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]';
        dcmHead.FrameOfReferenceUID = '';
        dcmHead.TemporalPositionIdentifier = i;
        dcmHead.NumberOfTemporalPositions = obj.nr_frames;
        dcmHead.TemporalResolution = obj.timeperframe;
        dcmHead.ImagesInAcquisition = obj.NO_SLICES;
        dcmHead.SliceLocation = (z-round(obj.NO_SLICES/2))*(obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE);
        dcmHead.ImageComments = '';
        dcmHead.TemporalPositionIndex = uint32([]);
        dcmHead.SamplesPerPixel = 1;
        dcmHead.PhotometricInterpretation = 'MONOCHROME2';
        dcmHead.PlanarConfiguration = 0;
        dcmHead.Rows = dimy;
        dcmHead.Columns = dimx;
        dcmHead.PixelSpacing = [pixely pixelx]';
        dcmHead.PixelAspectRatio = obj.aspectratio;
        dcmHead.BitsAllocated = 16;
        dcmHead.BitsStored = 15;
        dcmHead.HighBit = 14;
        dcmHead.PixelRepresentation = 0;
        dcmHead.PixelPaddingValue = 0;
        dcmHead.RescaleIntercept = 0;
        dcmHead.RescaleSlope = 1;
        dcmHead.HeartRate = 0;
        dcmHead.NumberOfSlices = obj.NO_SLICES;
        dcmHead.CardiacNumberOfImages = 1;
        dcmHead.MRAcquisitionType = convertStringsToChars(obj.dataType);
        dcmHead.ScanOptions = 'CG';
        dcmHead.BodyPartExamined = '';

        dicomHeader = dcmHead;

    end % Generate dicom header


    function dicomHeader = generateDicomheaderMRDflow(fname,fileCounter,j,k,z)

        % GENERATES DICOM HEADER FOR EXPORT OF FLOW IMAGES

        acq_dur = obj.nr_frames * obj.timeperframe;   % acquisition time in seconds

        if obj.PHASE_ORIENTATION == 1
            pixelx = obj.aspectratio*obj.FOV/dimx;
            pixely = obj.FOV/dimy;
        else
            pixelx = obj.FOV/dimx;
            pixely = obj.aspectratio*obj.FOV/dimy;
        end

        dt = datetime(obj.date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
        year = num2str(dt.Year);
        month = ['0',num2str(dt.Month)]; month = month(end-1:end);
        day = ['0',num2str(dt.Day)]; day = day(end-1:end);
        date = [year,month,day];

        hour = ['0',num2str(dt.Hour)]; hour = hour(end-1:end);
        minute = ['0',num2str(dt.Minute)]; minute = minute(end-1:end);
        seconds = ['0',num2str(dt.Second)]; seconds = seconds(end-1:end);
        time = [hour,minute,seconds];

        dcmHead.Filename = fname;
        dcmHead.FileModDate = obj.date;
        dcmHead.FileSize = dimy*dimz*2;
        dcmHead.Format = 'DICOM';
        dcmHead.FormatVersion = 3;
        dcmHead.Width = dimx;
        dcmHead.Height = dimy;
        dcmHead.BitDepth = 15;
        dcmHead.ColorType = 'grayscale';
        dcmHead.FileMetaInformationGroupLength = 178;
        dcmHead.FileMetaInformationVersion = uint8([0, 1])';
        dcmHead.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.TransferSyntaxUID = '1.2.840.10008.1.2.1';
        dcmHead.ImplementationClassUID = '1.2.826.0.9717382.3.0.3.6.0';
        dcmHead.ImplementationVersionName = 'OFFIS_DCMTK_360';
        dcmHead.SpecificCharacterSet = 'ISO_IR 100';
        dcmHead.ImageType = 'DERIVED\4D-MRI\';
        dcmHead.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.StudyDate = date;
        dcmHead.SeriesDate = date;
        dcmHead.AcquisitionDate = date;
        dcmHead.StudyTime = time;
        dcmHead.SeriesTime = time;
        dcmHead.AcquisitionTime = 1000*(k-1)*obj.timeperframe;
        dcmHead.ContentTime = time;
        dcmHead.Modality = 'MR';
        dcmHead.Manufacturer = 'MR Solutions Ltd';
        dcmHead.InstitutionName = 'Amsterdam UMC';
        dcmHead.InstitutionAddress = 'Amsterdam, Netherlands';
        dcmHead.ReferringPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.ReferringPhysicianName.GivenName = '';
        dcmHead.ReferringPhysicianName.MiddleName = '';
        dcmHead.ReferringPhysicianName.NamePrefix = '';
        dcmHead.ReferringPhysicianName.NameSuffix = '';
        dcmHead.StationName = 'MRI Scanner';
        dcmHead.StudyDescription = 'XD-data';
        dcmHead.SeriesDescription = '';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.GivenName = '';
        dcmHead.PhysicianOfRecord.MiddleName = '';
        dcmHead.PhysicianOfRecord.NamePrefix = '';
        dcmHead.PhysicianOfRecord.NameSuffix = '';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PerformingPhysicianName.GivenName = '';
        dcmHead.PerformingPhysicianName.MiddleName = '';
        dcmHead.PerformingPhysicianName.NamePrefix = '';
        dcmHead.PerformingPhysicianName.NameSuffix = '';
        dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianReadingStudy.GivenName = '';
        dcmHead.PhysicianReadingStudy.MiddleName = '';
        dcmHead.PhysicianReadingStudy.NamePrefix = '';
        dcmHead.PhysicianReadingStudy.NameSuffix = '';
        dcmHead.OperatorName.FamilyName = 'manager';
        dcmHead.AdmittingDiagnosesDescription = '';
        dcmHead.ManufacturerModelName = 'MRS7024';
        dcmHead.ReferencedSOPClassUID = '';
        dcmHead.ReferencedSOPInstanceUID = '';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.DerivationDescription = '';
        dcmHead.FrameType = '';
        dcmHead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PatientID = '01';
        dcmHead.PatientBirthDate = date;
        dcmHead.PatientBirthTime = '';
        dcmHead.PatientSex = 'F';
        dcmHead.OtherPatientID = '';
        dcmHead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.OtherPatientName.GivenName = '';
        dcmHead.OtherPatientName.MiddleName = '';
        dcmHead.OtherPatientName.NamePrefix = '';
        dcmHead.OtherPatientName.NameSuffix = '';
        dcmHead.PatientAge = '1';
        dcmHead.PatientSize = [];
        dcmHead.PatientWeight = 0.0300;
        dcmHead.Occupation = '';
        dcmHead.AdditionalPatientHistory = '';
        dcmHead.PatientComments = '';
        dcmHead.BodyPartExamined = '';
        dcmHead.SequenceName = obj.PPL;
        dcmHead.SliceThickness = obj.SLICE_THICKNESS;
        dcmHead.KVP = 0;
        dcmHead.RepetitionTime = obj.TR;
        dcmHead.EchoTime = obj.TE;                 % ECHO TIME
        dcmHead.InversionTime = 0;
        dcmHead.NumberOfAverages = obj.NO_AVERAGES;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.MagneticFieldStrength = obj.field_strength;
        dcmHead.SpacingBetweenSlices = obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE;
        dcmHead.EchoTrainLength = obj.NO_ECHOES;
        dcmHead.DeviceSerialNumber = '0034';
        dcmHead.PlateID = '';
        dcmHead.SoftwareVersion = '1.0.0.0';
        dcmHead.ProtocolName = '';
        dcmHead.SpatialResolution = [];
        dcmHead.TriggerTime = 1000*(k-1)*obj.timeperframe;  % trigger time in ms
        dcmHead.DistanceSourceToDetector = [];
        dcmHead.DistanceSourceToPatient = [];
        dcmHead.FieldofViewDimensions = [obj.FOV obj.aspectratio*obj.FOV obj.SLICE_THICKNESS];
        dcmHead.ExposureTime = [];
        dcmHead.XrayTubeCurrent = [];
        dcmHead.Exposure = [];
        dcmHead.ExposureInuAs = [];
        dcmHead.FilterType = '';
        dcmHead.GeneratorPower = [];
        dcmHead.CollimatorGridName = '';
        dcmHead.FocalSpot = [];
        dcmHead.DateOfLastCalibration = '';
        dcmHead.TimeOfLastCalibration = '';
        dcmHead.PlateType = '';
        dcmHead.PhosphorType = '';
        dcmHead.AcquisitionMatrix = uint16([dimx 0 0 dimy])';
        dcmHead.FlipAngle = obj.flipAngleArray(j);           % FLIP ANGLES
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.CassetteOrientation = 'PORTRAIT';
        dcmHead.CassetteSize = '25CMX25CM';
        dcmHead.ExposuresOnPlate = 0;
        dcmHead.RelativeXrayExposure = [];
        dcmHead.AcquisitionComments = '';
        dcmHead.PatientPosition = 'HFS';
        dcmHead.Sensitivity = [];
        dcmHead.FieldOfViewOrigin = [];
        dcmHead.FieldOfViewRotation = [];
        dcmHead.AcquisitionDuration = acq_dur;
        dcmHead.StudyInstanceUID = dcmid(1:18);
        dcmHead.SeriesInstanceUID = [dcmid(1:18),'.',num2str(obj.filename)];
        dcmHead.StudyID = '01';
        dcmHead.SeriesNumber = obj.filename;
        dcmHead.AcquisitionNumber = 1;
        dcmHead.InstanceNumber = fileCounter;
        dcmHead.ImagePositionPatient = [-obj.FOV/2 -(obj.aspectratio*obj.FOV/2) (z-round(obj.NO_SLICES/2))*(obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE)]';
        dcmHead.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]';
        dcmHead.FrameOfReferenceUID = '';
        dcmHead.TemporalPositionIdentifier = k;
        dcmHead.NumberOfTemporalPositions = obj.nr_frames;
        dcmHead.TemporalResolution = obj.timeperframe;
        dcmHead.ImagesInAcquisition = obj.NO_SLICES;
        dcmHead.SliceLocation = (z-round(obj.NO_SLICES/2))*(obj.SLICE_SEPARATION/obj.SLICE_INTERLEAVE);
        dcmHead.ImageComments = '';
        dcmHead.TemporalPositionIndex = uint32([]);
        dcmHead.SamplesPerPixel = 1;
        dcmHead.PhotometricInterpretation = 'MONOCHROME2';
        dcmHead.PlanarConfiguration = 0;
        dcmHead.Rows = dimy;
        dcmHead.Columns = dimx;
        dcmHead.PixelSpacing = [pixely pixelx]';
        dcmHead.PixelAspectRatio = obj.aspectratio;
        dcmHead.BitsAllocated = 16;
        dcmHead.BitsStored = 15;
        dcmHead.HighBit = 14;
        dcmHead.PixelRepresentation = 0;
        dcmHead.PixelPaddingValue = 0;
        dcmHead.RescaleIntercept = 0;
        dcmHead.RescaleSlope = 1;
        dcmHead.HeartRate = 0;
        dcmHead.NumberOfSlices = obj.NO_SLICES;
        dcmHead.CardiacNumberOfImages = 1;
        dcmHead.MRAcquisitionType = convertStringsToChars(obj.dataType);
        dcmHead.ScanOptions = 'CG';
        dcmHead.BodyPartExamined = '';

        dicomHeader = dcmHead;

    end % Generate dicom header flow



% ---------- END OF FUNCTIONS ----------------

end % exportDicomMRD