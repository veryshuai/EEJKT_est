Cluster=parcluster('LionX');
% Cluster.ResourceTemplate='-l nodes=2:ppn=3 -l walltime=100:00:00 -l pmem=1gb -q lionxg-econ';
% batch(Cluster,'batch_run','Matlabpool',5);

Cluster.ResourceTemplate='-l nodes=7:ppn=3 -l walltime=100:00:00 -l pmem=1gb -q lionxg-econ';
batch(Cluster,'batch_run','Matlabpool',20);