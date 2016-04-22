function [Tv Nv] = torus_tansp(xin)
%We build the tangent and normal spaces to the amplitude torus, for a 
%n-dimensional array,  at xin, which represents the point in image space
%
%This is done in two steps, the first step is to
%compute the list of indices and phases for the eventual tangent and normal
%vectors.
%The vectors themselves are built as a second step.
%
%Compute the Fourier representation
xah=fftn(xin);
fsz=sizer(xah); %The reduced size.
%
fsz0=size(xah); %The usual size.
%
%We determine the dimension, n of the "ambient" space for the data.
%
n=length(fsz);
fsl=prod(fsz(1:n)); %The number of elements in each FFT
M=sqrt(fsl); %Normalizing factor for the inverse FFT
%
%We get pointers into the flattened indices needed to build the Fourier
%representation of the tangent and normal spaces
[TInd NInd NInd0] = torus_tansp_ind(xah);
%
%Initialize an template for  basis vector of the normal and tangent and normal spaces
%in the Fourier representation.
%

%

%This gives and o/n basis for the real tangent space at xa  
%on the torus defined by the magnitude data of xah
tsz=size(TInd);
nsz=size(NInd);
n0=length(NInd0);
%
%We build the flattened template for the Fourier tangent/normal spaces
hFv=zeros(fsl,1);
%
for j=1:tsz(1)
%     tic
    hFv(TInd{j,1})=1i*TInd{j,3};
    hFv(TInd{j,2})=-1i*conj(TInd{j,3});
    Fv=M*ifftn(reshape(hFv,fsz0));%Maybe this should be done explicitly as there 
                          %are only 2 non-zero Fourier coefficients.
    Tv(1:fsl,j)=Fv(:); 
    %zero out locations for next iterate:                      
    hFv(TInd{j,1})=0;
    hFv(TInd{j,2})=0;   
%    if (mod(j,128) == 0)
%        j
%        length(Tv(:))
%        toc
%    end
end
%Reshape data into a sequence of n-diml arrays
Tv=reshape(Tv,[fsz,tsz(1)]);
%We build the special normal space first
%
for j=1:n0
    hFv(NInd0(j))=1;    %These are directions where the FT is real.
    Fv=M*ifftn(reshape(hFv,fsz0));%Maybe this should be done explicitly as there 
                                 %is only 1 non-zero Fourier coefficient.
    Nv(1:fsl,j)=Fv(:); 
    %zero out location for next iterate: 
    hFv(NInd0(j))=0;
end
%We add in the rest of the normal space
for j=1:nsz(1)
    hFv(NInd{j,1})=NInd{j,3};
    hFv(NInd{j,2})=conj(NInd{j,3});
    Fv=M*ifftn(reshape(hFv,fsz0));
    Nv(1:fsl,j+n0)=Fv(:);       %Maybe this should be done explicitly as there 
                                %are only 2 non-zero Fourier coefficients.
    %zero out locations for next iterate: 
    hFv(NInd{j,1})=0;
    hFv(NInd{j,2})=0;
end
%Reshape the data into a sequence of n-diml arrays
Nv=reshape(Nv,[fsz,nsz(1)+n0]);
end

function [TInd NInd NInd0] = torus_tansp_ind(xah)
%xah is a point on the amplitutde torus in the Fourier representation.
%The function outputs a list of the flattened indices where non-zero entries go to
%define orthonormal bases for the normal and tangent spaces to the
%amplitude torus at the base point xah:
%
%TInd{i,j} is a cell structure with i=1:3 and j indexes the tangent
%directions to the amplitude torus
%TInd{1,j}, TInd{2,j} are flat indices in the array defining a tangent vector 
%where the phase is TInd{3,j} for 1 and its conjugate for 2
%NInd{i,j} are similarly defined.
%NInd0(j) are indices of special normal directions that arise due to
%symmetries forcing certain Fourier coefficients to be real. For example,
%the zero Fourier coefficient with index (1,1,...,1) is always real and
%defines a special normal direction.
%
sqrt2=sqrt(2);
%
szr=sizer(xah);  %The lengths along various dimensions
%
n=length(szr);    %The number of dimensions 
%
tol = 10^(-14); %If and entry of xah is less than tol*fmx we will declare 
                %it to be 0, and hence not define a tangent direction.
fmx=max(abs(xah(:))); % This is the maximum modulus of the Fourier coefficients
%
%We loop through an index into the flattened array of the same size as xah,
%building the tangent and normal spaces to the torus at the point xah.
%
tcnt = 0;   %Count the number of tangent dimensions
ncnt = 0;   %Count the number of normal dimensions
ncnt0 = 0;  %Count the number of special normal dimensions
%
%
lx=length(xah(:)); %The length of the flattened array
for j=1:lx
    %Find the conjugate index to i2s(j)=M+2-i2s(j) mod M
    cj=conjind(szr,j);
    if cj == j
        %These entries are real and so there is only a single
        %special normal direction
        ncnt0=ncnt0+1;
        NInd0(ncnt0)=j;
        if tcnt+ncnt+ncnt0 == lx  %Check to see if we're done, should not be needed
            break
        end
    else
        %Now check to see if we've already found the pair (j,cj): if j and cj
        %are not equal then we can just keep the pairs where j < cj.
        if cj<j
            continue
        end
        %if cj =/= j or cj < j then we get both a tangent and normal direction
         xph=xah(j); %Find the Fourier component with this flat index
        if abs(xph)<tol*fmx %If an entry is too small then don't create
                            %a tangent  direction, but there are 2 normal
                            %directions
           ncnt=ncnt+1;
           NInd{ncnt,1}=j;
           NInd{ncnt,2}=cj;
           NInd{ncnt,3}=1;
           ncnt=ncnt+1;
           NInd{ncnt,1}=j;
           NInd{ncnt,2}=cj;
           NInd{ncnt,3}=1i;
           if tcnt+ncnt+ncnt0 >= lx  %Check to see if we're done, should not be needed
                 break
           end
           %
        else            %The Fourier component is non-zero and we create 
                        %a tangent and a normal direction.
           ncnt=ncnt+1;
           NInd{ncnt,1}=j;
           NInd{ncnt,2}=cj;
           NInd{ncnt,3}=xph/(sqrt2*abs(xph)); %Record phase
           tcnt=tcnt+1;
           TInd{tcnt,1}=j;
           TInd{tcnt,2}=cj;
           TInd{tcnt,3}=xph/(sqrt2*abs(xph)); %Record phase
           if tcnt+ncnt+ncnt0 >= lx  %Check to see if we're done, should not be needed
               break
           end
        end
    end
end
end

function sz=sizer(A)
%This returns the reduced size of an array, if the dimension is at least 2
%otherwise it returns the length.
sz0=size(A);
if (length(sz0) == 2) & (min(sz0) == 1)
    sz=length(A);
else
    sz=sz0;
end
end


    
