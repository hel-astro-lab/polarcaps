# #!/bin/zsh

declare -a conf_arr=(
#"turb_test.ini"
#"x8s10r3.ini"
#"x8s10r3a1.ini"
## v1
#"x16m30_p8np4c1_s100d03_r10_a05w08g06/x16s10.ini"
#"x16m30_p8np4c1_s10d03_r10_a05w08g06/x16s10.ini"
#"x16m30_p8np4c1_s10d03_r3_a05w08g06/x16s10.ini"
#"x16m30_p8np4c1_s10d03_r3_a08w08g06/x16s10.ini"
#"x16m30_p8np4c1_s10d03_r10_a08w08g06/x16s10.ini"
#"x16m30_p8np4c1_s10d03_r3_a1w08g06/x16s10.ini"
#"x16m30_p8np8c1_s10d03_r3_a05w08g06/x16s10.ini"
#"x16s10r3a1/x16s10r3a1.ini"
#v2
#"x16m30_p8np4c1_s100d03_r10_a05w08g06_v2/x16s10.ini"
#"x16m30_p8np4c1_s10d03_r10_a07w08g06_v2/x16s10.ini"
#"x32m30_p8np4c1_s10d03_r10_a07w08g06_v2/x32.ini"
#"x32m30_p8np4c1_s10d03_r10_a07w08g06_v3/x32.ini"
#"x32m40_p8np3c1_s10d03_r20_a07w08g06_v3/x32.ini"
#"x32m40_p6np4c1_s10d03_r10_a07w08g06_v4/x32.ini"
#
#"x16m40_p8np0c1_s10d03_r10_a07w08g06_v6/x16.ini"
#"x16m40_p8np2c1_s10d03_r10_a07w08g06_v6/x16.ini"
#"x16m40_p8np4c1_s10d03_r10_a07w08g06_v6/x16.ini"
#"x16m40_p8np6c1_s10d03_r10_a07w08g06_v6/x16.ini"
#"x16m40_p8np8c1_s10d03_r10_a07w08g06_v6/x16.ini"
#"x16m40_p32np0c1_s10d03_r10_a07w08g06_v6/x16.ini"
"cygx1-hard-v23-currleak2/cyg_hard.ini"
)




declare -a single_scripts_arr=(
"plot_ene_evolution.py"
"plot_prtcl_spec.py"
)


declare -a lap_scripts_arr=(
"plot3d.py --var jz" 
"plot3d.py --var bperp" 
"plot3d.py --var rho" 
#"plot3d.py --var rho" 
#"plot3d.py --var jz" 
#"plot3d.py --var bperp" 
#"plot3d.py --var bz" 
#"plot_fld_2d_panel.py --var rho" 
#"plot_fld_2d_panel.py --var jx" 
#"plot_fld_2d_panel.py --var jy" 
#"plot_fld_2d_panel.py --var ey" 
#"plot3d.py --var mpi" 
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
for ((lap=0;lap<=10000;lap+=10))
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

