#!/bin/bash
myArray=("2b_extract_sign.R" "3a_hyprcoloc_gtex.R" "3b_computehyprcoloc_RNAProtBMI.R" "3c_vep.R" "3d_ztest_ukb_giant.R" "4e_forarnaud.R" "5_tissular_specificity.R" "5a_manhattan_plot.R" "5c_hyprcoloc_plot.R" "5e_Venndiagramm.R" "5f_magma_histogramm.R" "5f_PhenoGram.R"  "5g_create_tables.R")                      

for str in ${myArray[@]}; do
  chmod u+x ./$str 
done

echo "Initializing 2b_extract_sign" && ./2b_extract_sign.R &&
 echo "Initializing 3a_hyprcoloc_gtex.R" && ./3a_hyprcoloc_gtex.R && 
 echo "Initializing 3b_computehyprcoloc_RNAProtBMI.R" && ./3b_computehyprcoloc_RNAProtBMI.R &&
 echo "Initializing 3c_vep.R" && ./3c_vep.R &&
 echo "Initializing 3d_ztest_ukb_giant.R" && ./3d_ztest_ukb_giant.R && 
 echo "Initializing 4e_forarnaud.R" && ./4e_forarnaud.R && 
 echo "Initializing 5_tissular_specificity.R" && ./5_tissular_specificity.R &&
 echo "Initializing 5a_manhattan_plot.R" && ./5a_manhattan_plot.R && 
 echo "Initializing 5c_hyprcoloc_plot.R" && ./5c_hyprcoloc_plot.R &&
echo "Initializing 5e_Venndiagramm.R" && ./5e_Venndiagramm.R &&
echo "Initializing 5f_magma_histogramm.R" && ./5f_magma_histogramm.R &&
echo "Initializing 5f_PhenoGram.R" && ./5f_PhenoGram.R &&
echo "Initializing 5g_create_tables.R" && ./5g_create_tables.R &&
echo "The master script finished without errors"
