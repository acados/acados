#!/bin/bash
nfunctions=0 
find ./build/acados/CMakeFiles/acados.dir -type d | while read directory 
do  
    printf "\n\n----------------------------------------------------------\n"
    printf "directory: $directory"
    printf "\n----------------------------------------------------------\n"
    for file in $directory/* ; do
        if [ ${file: -2} == ".o" ]
        then
            printf "\n"
            echo "-> file: $file"
            printf "\n"
            nm ${file} | grep "T " 
            printf "\n"
            echo "number of functions to be documented: "
            nm ${file} | grep "T " | grep -v " _" | wc -l
            read x words chars filename <<< $(nm ${file} | grep "T " | grep -v " _" | wc )
            nfunctions=$((nfunctions+x))
            printf "\n"
        fi
    done <<< $(find tmp -type f)
    printf "\n\ntotal number of functions: $nfunctions\n\n"
done  
