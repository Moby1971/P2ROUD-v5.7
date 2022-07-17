function  res = Emat_zyxt(mask,b1)

%   res = Emat_zyxt(mask,b1)
%
%
%	Implementation of parallel MRI encoding matrix for dynamic MRI data
%	
%	input:
%			mask : ky-kx-t sampling mask (Nx,Ny,Nz,Nt)
%           b1 : coil sensitivity maps (Nx,Ny,Nz,Nc)
%
%	output: the operator
%

res.adjoint = 0;
res.mask = mask;
res.b1 = b1;
res = class(res,'Emat_zyxt');