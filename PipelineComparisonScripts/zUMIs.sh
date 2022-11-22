#!/bin/bash/

timefile="/media/draco/lucask/AlgorithmTest/zUMI/timelog.txt"

echo -e "START: $(date)\n" >> ${timefile}

conda activate zUMI
echo $CONDA_DEFAULT_ENV

/media/draco/lucask/zUMIs/zUMIs.sh -c -y /media/draco/lucask/AlgorithmTest/zUMI/zUMI_benchmark_test.yaml

conda deactivate
echo $CONDA_DEFAULT_ENV

echo -e "STOP: $(date)\n" >> ${timefile}

