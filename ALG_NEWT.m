function [D,recon,info]=ALG_NEWT(D,u,algopts)
%This function uses linear programing/or least squares
%to iteratively improve the guess, D.x for
%an image
% D is the updated guess
% info.Tv is an orthonormal basis for the tangent space to the magnitude torus
% info.resid is the residual between D.x_in and D.x_out
%The nearest point on the magnitude torus (in Fourier space) is found here:
%
%For calls to pi_A
opts.u=u;
%Remember the input image for computing the residual
xin=D.x;
%Find the nearest point to xin on the magnitude torus
xa=pi_A(D.x,opts);
%
%Compute the tangent space to the torus at the xa
%
[Tv Nv] = torus_tansp(xa);
%
%Return the o/n basis for the tangent space for latter analysis
info.Tv=Tv;
%Tv and Nv are cell arrays of vectors that need to be converted to ordinary
%arrays
%We want to know if there is a vector w so that xa+Tv*w has positive
%entries. This can be rephrased as a linear program: we look at vectors v
%such that Nv'*v=Nv'*xa, with non-negative coefficients and we try to find
%the one with the smallest L1 norm in the complement of the known region 
%where the correction solution is known to be supported.
%
%We flatten the support mask
fs0=algopts.support_mask(:);
%We flatten the data
fxa=xa(:);
%and the representation of the normal space to Tv_xa
szn=size(Nv);
nn=length(szn); %Dimension of the ambient space + 1
mm=prod(szn(1:(nn-1)));
%Flattened representation of the o/n basis for the normal space
fNv=reshape(Nv,[mm,szn(nn)]);
%
%These equality constraints define the tangent space to A at xa
Aeq=fNv';
beq=fNv'*fxa;
%
%If we want a non-negative solution
if algopts.positivity == 1
    lb=zeros(size(fxa));
    %More or less arbitrary upper bound:
    ub=4*abs(sum(fxa))*ones(size(fxa));
else %No upper and lower bounds
    lb=-ones(size(fxa))*inf;
    ub=ones(size(fxa))*inf;
end

%We want to minimize the square sum of the coordinates over the complement
%of the KNOWN support, and we may also include a penalty term
%of size eps_reg to prevent wandering in the support subspace.
eps_reg=0;
if algopts.support == 1
    eps_reg=algopts.eps_reg;
end
    Sc=eye(length(fs0))-diag(fs0)+eps_reg*diag(fs0);
    vc=zeros(size(fs0));

%Here we use the linear least squares solver to impose the positivity and/or
%support conditions
[xout,fval]=lsqlin(Sc,vc,[],[],Aeq,beq,lb,ub);
if nn > 2
    D.x=reshape(xout,[szn(1:(nn-1))]);
else
    D.x=xout;
end
recon=pi_A(D.x,opts);
info.resid=norm(xin-D.x)/norm(xin);
end

