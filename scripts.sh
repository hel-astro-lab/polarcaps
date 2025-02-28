# #!/bin/zsh

declare -a conf_arr=(
#"test.ini"
"test1d.ini"
)


declare -a single_scripts_arr=(
#"plot_ene_evolution.py"
#"plot_prtcl_spec.py"
)


declare -a lap_scripts_arr=(
#"plot3d_pyvista.py --var bvec --view side" 
#"plot3d_pyvista.py --var bvec --view tilt" 
#"fig_cascade.py "
#"fig_gap_TS_RCF.py"
"fig_gap.py"
#"fig_cascade.py"
#"plot3d_pyvista.py --var bvec --view top" 
#"plot_pulsar_lc.py "
)




#-------------------------------------------------- 
# Single lap scripts
for c in "${conf_arr[@]}"
do
    echo "conf filename is: $c"
    for s in "${single_scripts_arr[@]}"
    do
        echo "script is: $s"
        python3 $s --conf $c
    done
done

#-------------------------------------------------- 
# multilap scripts

#for lap in {0..200..10}
for ((lap=0;lap<=30000;lap+=200))
do
    echo "lap is $lap"
    for c in "${conf_arr[@]}"
    do
        echo "conf filename is: $c"
        for s in "${lap_scripts_arr[@]}"
        do
            echo "script is: $s"
            python3 $s --conf $c --lap $lap
        done
    done
done

