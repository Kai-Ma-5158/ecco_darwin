#!/bin/bash
#SBATCH --job-name=LLC90-new         # -J, 作业名称，方便yhq时查看区分作业信息
#SBATCH --partition=cp6           # -p, 资源分区名称, 在HPC2系统上为TH_HPC2
#SBATCH --nodes=2                # -N, 指定作业使用的节点数
#SBATCH --ntasks=96              # -n, 提交作业的总进程数

yhrun -N 2 -n 96 -p cp6 ./mitgcmuv 

#yhrun -N 15 -n 767 -p cp6 ./mitgcmuv > LLC540_log2.in

#if use debug
#yhrun -N 1 -n 12 -p debug2 ./romsM ./ocean_pre.in > pre_log.in
