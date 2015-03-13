% This file runs batch_run.m in batch mode. Number of nodes should be 

% Matlabpool plus 1.



Cluster=parcluster('LionX');

Cluster.ResourceTemplate='-V -l nodes=26 -l walltime=90:00:00 -l pmem=10gb -q lionxg-econ';

job = batch(Cluster,'batch_run','Matlabpool',26,'CaptureDiary',true);




