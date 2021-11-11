#!/bin/bash

CPU_SOCKET=`lscpu|awk '/Socket\(s\):/{print $2}'`
NUMA_NODES=`lscpu|awk '/NUMA node\(s\):/{print $3}'`
SYSCORES=`lscpu|awk '{if($1=="CPU(s):") {print $2}}'`

echo "#define TPSM_CPUSOCKET  ($CPU_SOCKET)"
echo "#define TPSM_NUMANODES  ($NUMA_NODES)"
echo "#define TPSM_SYSCORES   ($SYSCORES)"
