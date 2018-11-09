find ./acados -type f -not -path "*Makefile*" -exec sed -i 's/	/    /g' {} \;
find ./interfaces/acados_c -type f -not -path "*Makefile*" -exec sed -i 's/	/    /g' {} \;
