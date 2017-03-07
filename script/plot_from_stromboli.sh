#!/bin/bash

echo "plot_from_stromboli.sh starts with argument: " $@

if [ $# = 1 ]; then
    JSONFiles=$@
    postfix=
else
    postfix=${@:$#}
    JSONFiles=${*%${!#}}
fi

echo "-- JSONFiles list is: " $JSONFiles
echo "-- postfix list is: " $postfix

echo "-- Loading python 2.7"
source activate py27

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

echo "-- Looping over listed input JSONFiles"
for JSONFile in $JSONFiles
do
    prefix=$(python parse_json.py $JSONFile prefix)

    echo "-- Copying datas to local directory"
    scpf stromboli home/tauzin/Workspace/lbm_solver/output/outputPy/${prefix}_* ../output/outputPy/
    scpf stromboli home/tauzin/Workspace/lbm_solver/output/outputPy/${prefix}/ ../output/outputPy/

    echo "-- Copying plots to local directory"
    scpf stromboli home/tauzin/Workspace/lbm_solver/output/outputPlot/${prefix}_* ../output/outputPlot/
    scpf stromboli home/tauzin/Workspace/lbm_solver/output/outputAnimation/${prefix}_* ../output/outputAnimation/
    scpf stromboli home/tauzin/Workspace/lbm_solver/output/outputPlot/${prefix}/ ../output/outputPlot/
    scpf stromboli home/tauzin/Workspace/lbm_solver/output/outputAnimation/${prefix}/ ../output/outputAnimation/


done

screen -dm bash -c 'cd Workspace/lbm_solver/script && ./plot.sh ${JSONFiles} ${postfix} &> plot_${postfix}.out'"

echo "-- Opening generated plots in nautilus"
nautilus ../output/outputPlot
nautilus ../output/outputAnimation

echo "lbm_on_stromboli.sh ends"
