function [D,recon,info]=ALG_AB_HIO(D,u,algopts)

algopts.u=u;
f0=D.x;
pi_B_f0=pi_B(f0,algopts);
proj1=pi_A(2*pi_B_f0-f0,algopts);
proj2=pi_B_f0;
f1=f0+proj1-proj2;
D.x=f1;

recon=pi_B(D.x,algopts);
info.resid=compute_residual(proj1,proj2);

end