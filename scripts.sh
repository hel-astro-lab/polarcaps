# #!/bin/zsh

declare -a conf_arr=(
"test3d.ini"
)




declare -a single_scripts_arr=(
#"plot_ene_evolution.py"
#"plot_prtcl_spec.py"
)


declare -a lap_scripts_arr=(
"plot3d_pyvista.py --var bvec --view side" 
#"plot3d_pyvista.py --var bvec --view top" 
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
for ((lap=0;lap<=2000;lap+=20))
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

