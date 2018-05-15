#!/usr/bin/env python
from subprocess import run, PIPE, STDOUT
from pathlib import Path
from sys import argv, path

cluster = argv[1]
scaling = argv[2]
num_simulations = int(argv[3])
domain_len = int(argv[4])

# write_step = argv[4] # for persistant kernels

# Iterates over processes with powers of 2
for procs in [ 2**j for j in range(0, num_simulations) ]:
  x_len = domain_len
  y_len = domain_len
  z_len = domain_len
  if scaling == 'weak':
    x_len = procs * domain_len
  
  lbm_postfix = 'benchmark_{scaling}_{procs}'.format(**locals())
  scaling_cmd = './batch_{cluster}.sh {procs} {x_len} {y_len} {z_len} {lbm_postfix}'.format(**locals())
  run(strong_scaling_cmd, cwd=Path(path[0]) shell=True)
