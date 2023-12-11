% ---------------------------------------------------------------------------------
% Combine images from different coils
% ---------------------------------------------------------------------------------
function im2 = coilCombine(im1)

im1 = permute(im1,[1 2 3 5 4]);
[sx,sy,sz,~,N] = size(im1);

% Loop over slices
im2 = complex(zeros(sx,sy,sz,1,N));
for kz=1:sz
    im2(:,:,kz,1,:) = coilCombine(im1(:,:,kz,:,:));
end

im2 = permute(im2,[1 2 3 5 4]);

% Single slice coil combine
    function im2 = coilCombine(im1)

        im1 = permute(im1,[1 2 5 4 3]);

        % Get image dimensions and set filter size
        [sx,sy,N,C] = size(im1);
        filtsize = 7;

        % Initialize
        im2 = complex(zeros(sx,sy,1,1,N));
        Rs = complex(zeros(sx,sy,C,C));

        % Get correlation matrices
        for kc1=1:C
            for kc2=1:C
                for kn=1:N
                    Rs(:,:,kc1,kc2) = Rs(:,:,kc1,kc2) + filter2(ones(filtsize),im1(:,:,kn,kc1).*conj(im1(:,:,kn,kc2)),'same');
                end
            end
        end

        % Compute and apply filter at each voxel
        for kx=1:sx
            for ky=1:sy
                [U,~] = svd(squeeze(Rs(kx,ky,:,:)));
                myfilt = U(:,1);
                im2(kx,ky,1,1,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);
            end
        end

    end


    end % coilCombine
