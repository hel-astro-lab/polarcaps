# #!/bin/zsh

declare -a conf_arr=(
#"test.ini"
#"test1d.ini"
#"rp_x1024m20_p2np24_v3e-06inje01injx0_b002N3e12_v6-int200-debug/rp.ini"
#"rp_x1024m20_p2np24_v3e-06inje01injx0_b002N3e12_v6-jm-1/rp.ini"
#"rp_x1024m20_p2np24_v3e-06inje01injx0_b002N3e12_v6-jm2/rp.ini"
#"rp_x1024m20_p2np24_v3e-06inje001injx0_b002N3e12_v7-jm1.5-ninj0.01/rp.ini"
#"rp_large_x2048m40_p2np12_v3e-06inje01injx0_b002N6e11_v6/rp-large.ini"
#"rp_large_x2048m40_p2np12_v3e-06inje0001injx0_b002N6e11_v7-ninj-3/rp-large.ini"
#"rp_large_x2048m40_p2np12_v3e-06inje001injx0_b002N6e11_v8-ninj-2-nojrot/rp-large.ini"
#"rp_x1024m20_p2np24_v3e-06inje001injx0_b002N3e12_v8-jm1.5-ninj0.01-nojrot/rp.ini"
#"rp_x1024m20_p2np128_v3e-06inje001injx0_b002N3e12_v8-jm1.5-ninj0.01-nojrot-flt/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_b002N3e12_v9-jm1.5/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc3_b002N2e12_v10/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v11/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v12/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v12-max100/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v14-nc/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v14-cap100/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v14-cap1000/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v15-c100/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v15-c1000/rp.ini"
#"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v15-nocap/rp.ini"
"rp_x1024m20_p2np32_v3e-06inje001injx0_jm15rc10_b002N2e12_v16-c2/rp.ini"
)


declare -a single_scripts_arr=(
#"plot_ene_evolution.py"
#"plot_prtcl_spec.py"
"fig_gap_dynamics.py"
)


declare -a lap_scripts_arr=(
#"plot3d_pyvista.py --var bvec --view side" 
#"plot3d_pyvista.py --var bvec --view tilt" 
#"fig_cascade.py "
#"fig_gap_TS_RCF.py"
#"fig_gap.py"
"fig_gap2.py"
#"fig_kspec.py"
#"fig_prtcl_profiles.py"
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
#for ((lap=230000;lap<=1000000;lap+=1000))
#for ((lap=28000;lap<=100000;lap+=2000))
for ((lap=0;lap<=100000;lap+=2000))
do
    echo "lap is $lap"
    for c in "${conf_arr[@]}"
    do
        echo "conf filename is: $c"
        for s in "${lap_scripts_arr[@]}"
        do
            echo "script: $s --- lap: $lap --- conf: $c"
            python3 $s --conf $c --lap $lap
        done
    done
done

