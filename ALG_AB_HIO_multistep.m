function [D,recon,info]=ALG_AB_HIO_multistep(D,u,algopts)

if (length(find(algopts.restart_steps==D.itnum))>0)
    D.x=pi_B(D.x,algopts);
end;
[D,recon,info]=ALG_AB_HIO(D,u,algopts);

end