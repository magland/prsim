function [D,recon,info]=ALG_BA_HIO(D,u,algopts)

algopts.u=u;
f0=D.x;
pi_A_f0=pi_A(f0,algopts);
proj1=pi_B(2*pi_A_f0-f0,algopts);
proj2=pi_A_f0;
f1=f0+proj1-proj2;
D.x=f1;

recon=pi_A(D.x,algopts);
info.resid=compute_residual(proj1,proj2);

end