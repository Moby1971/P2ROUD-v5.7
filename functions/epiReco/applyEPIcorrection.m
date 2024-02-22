function kSpaceCorrected = applyEPIcorrection(kSpaceGhost, kCenter, pCenter, method)

% Apply a first order even/odd k-space line correction

[kappa, phi] = findEPIshift(kSpaceGhost, kCenter, pCenter, method);

kSpaceCorrected = firstOrderPhaseCorr(kSpaceGhost, kappa, phi);





% -----------------------
% find EPI shift
% -----------------------

    function [kappa, phi] = findEPIshift(kSpace, kCenter, pCenter, method)

        [dimx,dimy] = size(kSpace);

        % SVD Indices
        if strcmp(method,'svd')

            blockX = 3; % kernel size
            blockY = 3; % kernel size

            nBlock = blockX * blockY;

            % Precaculate the mapping from kspTmp to the Mat matrix
            idxMap = zeros((dimx-blockX+1) * (dimy-blockY+1), nBlock);

            % These matrices define the offsets for the kernel
            kernelOffsetX = repmat(0:blockX-1, [blockY 1])';
            kernelOffsetY = repmat(0:blockY-1, [blockX 1]);

            count = 0;
            for kx = 1:dimx-blockX+1
                for ky = 1:dimy-blockY+1
                    count = count + 1;
                    kernelIndices = sub2ind([dimx,dimy], kernelOffsetX+kx, kernelOffsetY+ky); % JAM
                    idxMap(count,:) = kernelIndices(:);
                end
            end
        else
            idxMap = [];
        end

        % Create an anonymous function that you can use to pass parameters to the
        % optimization function
        anonFun = @(x)calculateMetric(x, kSpace, method, idxMap);

        % Starting values
        x0 = zeros(1,2);
        x0(1) = kCenter;
        x0(2) = pCenter;

        % Specify fitting options
        options = optimset(...
            'Display', 'off', ...
            'TolFun', 1E-6, 'TolX', 1E-6, ...
            'MaxFunEvals', 2000);

        % Begin Fitting
        [x] = fminsearch(anonFun, x0, options); % JAM
    
        % Save outputs
        kappa = x(1);
        phi = x(2);

        return;

        % This calculates the metric that needs to be minimized
        % The optimization parameters kappa and phi are arrayed in the parameter
        % matrix "x"

        function metric = calculateMetric(x, ksp, method, idxMap)
            kappaTmp = x(1);
            phiTmp = x(2);

            % Apply these corrections
            kspMod = firstOrderPhaseCorr(ksp,kappaTmp,phiTmp);
            metric = getMetric(kspMod, method, idxMap); % JAM
        end


        function [metricOut] = getMetric(kSpace, method, idxMap)

            switch method

                case 'ent'
                    % Entropy method, from:
                    % Clare S. Iterative Nyquist ghost correction for single and
                    %   multishot EPI using an entropy measure. ISMRM, 2003. p. 1041.
                    imgTmp = fftshift(fft2(ifftshift(kSpace)));
                    metricOut = entropy(mat2gray(double(abs(imgTmp))));

                case 'entSmooth'
                    % Entropy, plus smoothing
                    imgTmp = fftshift(fft2(ifftshift(kSpace)));
                    imgSmooth = medfilt2(abs(imgTmp));
                    metricOut = entropy(mat2gray(double(abs(imgSmooth))));

                case 'svd'
                    % SVD method, from:
                    % Peterson E, Aksoy M, Maclaren J, Bammer R. Acquisition?free
                    %   Nyquist ghost correction for parallel imaging accelerated EPI.
                    %   ISMRM 2015. p. 75.
                    minZone = 2; % first index of min zone
                    Mat = kSpace(idxMap);
                    S = svd(Mat);
                    tail = S(minZone:end);
                    %         tail = S(minZone+1:end); % HACK
                    metricOut = sum(tail); % The integral of the "tail" of the SVD vector

                otherwise 
                    
                    % Ghost-Object method, from:
                    % McKay JA, Moeller S, Zhang L, Auerbach EJ, Nelson MT, Bolan PJ.
                    %   Nyquist ghost correction of breast diffusion weighted imaging
                    %   using referenceless methods. MRM 2019 Apr;81(4):2624-2631
                    imgTmp = fftshift(fft2(ifftshift(kSpace)));
                    shift = circshift(imgTmp, [0,size(imgTmp,2)/2]);

                    metGhOb = abs(imgTmp./shift); % signal/shifted signal
                    metGhOb = medfilt2(metGhOb);
                    mGh = metGhOb(:);
                    mGh(isnan(mGh)) = [];
                    metric = mean(mGh);
                    metricOut = 1./metric; % Invert for minizmation problem

            end

        end % getMetric

    end % simplexSearch




% -----------------------
% First order phase correction
% -----------------------

    function newKspace = firstOrderPhaseCorr(oldKspace, kappa, phi)

        newKspace = complex(zeros(size(oldKspace,1),size(oldKspace,2)));  

        if mod(size(oldKspace,2),2) % If odd
            newKspace = complex(zeros(size(oldKspace,1),size(oldKspace,2)+1));
        end

        newKspace(1:size(oldKspace,1),1:size(oldKspace,2)) = oldKspace;

        % Apply phase correction
        newKspace(:,1:2:end) = newKspace(:,1:2:end) .* exp(+1j*phi);
        newKspace(:,2:2:end) = newKspace(:,2:2:end) .* exp(-1j*phi);

        % FT to image space in readout direction, apply phase correction
        phaseCorrection = complex(zeros(size(newKspace)));
        [dimx, dimy] = size(newKspace);
        xunits = (1:dimx) - (dimx/2) - 0.5;

        % Specify the even and odd traces
        phaseCorrection(:,1) = exp(+1j * 2*pi*kappa * xunits./(dimx));
        phaseCorrection(:,2) = exp(-1j * 2*pi*kappa * xunits./(dimx));

        % Replicate them over the full matrix
        phaseCorrection = repmat(phaseCorrection(:,1:2), [1, dimy/2]);

        % Apply phase correction
        fttData = fftshift(fft(ifftshift(newKspace,1)),1);
        newKspace = ifftshift(ifft(ifftshift(fttData.* phaseCorrection,1)  ,[],1),1);

    end % firstOrderPhaseCorr



end % applyEPIcorrection

