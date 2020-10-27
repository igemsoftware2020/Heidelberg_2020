#!/bin/sh
hhblits -cpu 4 -i "$1"/"$2".seq -d /beegfs/work/ws/hd_vu199-hh_suite_arnoldt-0/databases/uniclust30_2018_08/uniclust30_2018_08 -oa3m ./output/trrosetta/"$2".a3m -n 2
python ./trRosetta/network/predict.py -m ./trRosetta/model2019_07 ./output/trrosetta/"$2".a3m ./output/trrosetta/distances
for i in $(seq 1 $3)
do
    python ./trRosetta/trRosetta.py ./output/trrosetta/distances.npz "$1"/"$2".fasta ./output/trrosetta/"$2"-"$i".pdb
done
