# Executables
nucmer='nucmer'
deltafilter='delta-filter'
mafft=''
partition=''
raxml='raxmlHPC-MPI-AVX'
parse_examl='parse-examl'
examl='examl-AVX'

# SLURM settings
SLURMpartition='savio'
SLURMaccount='co_rosalind'
SLURMqos='rosalind_savio_normal'
SLURMnodes=4 # max number of nodes to utilize
SLURMcpus=20 # number of CPUs per node
SLURMmem=3000 # memory (in MBi) per CPU (should be at least 2000)

# Modules
python='python/2.7.14' # needs biopython, numpy, ete3
gcc='gcc/6.3.0'
mpi='openmpi/2.0.2-gcc'
