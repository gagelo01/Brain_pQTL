#!/bin/bash
myArray=(5a_plot_all_figures.R 5f_PhenoGram.R 5g_create_tables.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing 5a_plot_all_figures.R' && ./5a_plot_all_figures.R &&
echo 'Initializing 5f_PhenoGram.R' && ./5f_PhenoGram.R &&
echo 'Initializing 5g_create_tables.R' && ./5g_create_tables.R && echo 'The master script finished without errors'
