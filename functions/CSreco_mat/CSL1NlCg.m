function x = CSL1NlCg(x0,param)

% res = CSL1NlCg(param)
%
% Reconstruction of undersampled k-space MRI data
%
% Given the acquisition model y = E*x, and the sparsifying transform W,
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda2 * TV(x)
%

% starting point
x = x0;

% line search parameters
maxlsiter = 20;
alpha = 0.01;
beta = 0.6;
t0 = 1;

gradToll = 1e-3 ; %#ok<NASGU>
param.l1Smooth = 1e-15;

% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% inner iterations
for k=1:param.nite

    % backtracking line-search
    f0 = objective(x,dx,0,param);
    t = t0;
    f1 = objective(x,dx,t,param);
    lsiter = 0;

    while (f1 > (f0 - alpha*t*abs(g0(:)'*dx(:)))) && (lsiter<maxlsiter)
        lsiter = lsiter + 1;
        t = t * beta;
        f1 = objective(x,dx,t,param);
    end

    % control the number of line searches by adapting the initial step search
    if lsiter > 2, t0 = t0 * beta; end
    if lsiter < 1, t0 = t0 / beta; end

    % update x
    x = (x + t*dx);

    % conjugate gradient calculation
    g1 = grad(x,param);
    bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    g0 = g1;
    dx = -g1+bk*dx;

end

end



function res = objective(x,dx,t,param)

% L2-norm
w = param.E*(x+t*dx) - param.y;
L2Obj = w(:)'*w(:);

% L1-norm part
if param.L1Weight
    w = param.W*(x+t*dx);
    L1Obj = sum((conj(w(:)).*w(:) + param.l1Smooth).^(1/2));
else
    L1Obj=0;
end

% Total variation
w = param.TV*(x+t*dx);
TVObj = sum((w(:).*conj(w(:))).^(0.5));

% Objective function
res = L2Obj + param.TVWeight*TVObj + param.L1Weight*L1Obj;

end



function g = grad(x,param)

% L2-norm
L2Grad = 2.*(param.E'*(param.E*x-param.y));

% L1-norm part
if param.L1Weight
    w = param.W * x;
    L1Grad = param.W'*(w.*(w.*conj(w) + param.l1Smooth).^(-0.5));
else
    L1Grad=0;
end

% Total variation
w = param.TV*x;
TVGrad = param.TV'*(w.*(w.*conj(w)).^(-0.5));

% Complete gradient
g = L2Grad + param.TVWeight*TVGrad + param.L1Weight*L1Grad ;

end
