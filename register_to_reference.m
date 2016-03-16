function A=register_to_reference(A,A_ref)

N1=size(A_ref,1);
N2=size(A_ref,2);
N3=size(A_ref,3);
M1=ceil((N1+1)/2);
M2=ceil((N2+1)/2);
M3=ceil((N3+1)/2);

tmp1=abs(ifft0(fft0(A).*conj(fft0(A_ref))));
tmp2=abs(ifft0(conj(fft0(A)).*conj(fft0(A_ref))));
[max1,ii1]=max(tmp1(:));
[max2,ii2]=max(tmp2(:));

if (max1>max2)
    ii=ii1;
else
    ii=ii2;
    A=ifft0(conj(fft0(A)));
end;

[i1,i2,i3]=ind2sub(size(A),ii);
A=circshift(A,-(i1-M1),1);
A=circshift(A,-(i2-M2),2);
A=circshift(A,-(i3-M3),3);

if (sum(A(:).*A_ref(:))<0)
    A=-A;
end;

dxrange=-1:0.2:1;
dyrange=-1:0.2:1;
dzrange=-1:0.2:1;
if (N3==1) dzrange=0; end;
best_resid=inf; best_A=A;
for dz=dzrange
for dy=dyrange
for dx=dxrange
    Atry=fracshift(A,dx,dy,dz);
    resid=Atry-A_ref; resid=sum(resid(:).^2);
    if (resid<best_resid)
        best_A=Atry;
        best_resid=resid;
    end;
end;
end;
end;

A=best_A;

% [GX,GY,GZ]=ndgrid((0:N1-1)-M1,(0:N2-1)-M2,(0:N3-1)-M3);
% A1_CX=sum(A_ref(:).*GX(:))/sum(A_ref(:));
% A1_CY=sum(A_ref(:).*GY(:))/sum(A_ref(:));
% A1_CZ=sum(A_ref(:).*GZ(:))/sum(A_ref(:));
% 
% A2_CX=sum(A(:).*GX(:))/sum(A(:));
% A2_CY=sum(A(:).*GY(:))/sum(A(:));
% A2_CZ=sum(A(:).*GZ(:))/sum(A(:));
% 
% DX=A2_CX-A1_CX;
% DY=A2_CY-A1_CY;
% DZ=A2_CZ-A1_CZ;
% 
% A2=real(ifft0(fft0(A).*exp(2*pi*i*(DX*GX/N1+DY*GY/N2+DZ*GZ/N3))));

end

function X=fracshift(X,dx,dy,dz)
Xhat=fftn(X);
[N1,N2,N3]=size(X);
%N=7:[0,1,2, 3,-3,-2,-1] or N=8:[0,1,2,3, -4,-3,-2,-1]
M1a=floor((N1-1)/2); M1b=floor(N1/2);
kxrange=[0:M1a,-M1b:-1];
M2a=floor((N2-1)/2); M2b=floor(N2/2);
kyrange=[0:M2a,-M2b:-1];
M3a=floor((N3-1)/2); M3b=floor(N3/2);
kzrange=[0:M3a,-M3b:-1];
[kx,ky,kz]=ndgrid(kxrange,kyrange,kzrange);
theta=(kx*dx/N1+ky*dy/N2+kz*dz/N3)*2*pi;
Xhat=Xhat.*exp(i*theta);
X=real(ifftn(Xhat));
end
