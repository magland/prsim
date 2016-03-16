function test_prsim

%close all;

simopts.num_trials=10;
simopts.num_restarts=1;
simopts.num_iterations=500;
simopts.seed=1;

algopts.positivity=1;
algopts.support=0;

phopts.N=128;
phopts.oversamp=2;
phopts.num_disks=50;
phopts.k=2;
PHANTOM=@PH_1;

results_AB=prsim(@ALG_AB_HIO,PHANTOM,algopts,phopts,simopts);
results_BA=prsim(@ALG_BA_HIO,PHANTOM,algopts,phopts,simopts);

%algopts.restart_steps=100:100:simopts.num_iterations;
%prsim(@ALG_AB_HIO_multistep,PHANTOM,algopts,phopts,simopts);

end