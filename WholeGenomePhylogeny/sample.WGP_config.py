# Executables
nucmer='nucmer'
deltafilter='delta-filter'
mafft='mafft'
tiger='tiger2'
raxml='raxmlHPC-MPI-AVX'
parse_examl='parse-examl'
examl='examl-AVX'

# SLURM settings
SLURMpartition='compute'
SLURMaccount=None
SLURMqos=None
SLURMmaxtime=48 # time in hours
SLURMnodes=4 # max number of nodes to utilize
SLURMcpus=24 # number of CPUs per node
SLURMmem=3000 # memory (in MBi) per CPU (should be at least 2000)

# Large memory SLURM settings
LARGEpartition='large-shared'
LARGEaccount=None
LARGEqos=None
LARGEmaxtime=48 # time in hours
LARGEnodes=1 # max number of nodes to utilize
LARGEcpus=2 # number of CPUs per node
LARGEmem=48000 # memory (in MB) per CPU (should be at least 8000)


# Modules
modules = []
# needs biopython, numpy, ete3
# using personal anaconda python for this

