#!/bin/bash

source ~/.bashrc

set -e

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <pipeline_file> [rscript_command]"
    echo ''
    echo "The [rscript_command] parameter can be utilized to specify an architecture and a path to Rscript or a path to Rscript."
    echo "Example [rscript_command] parameter: \"arch -x86_64 /Users/michaelbertagna/anaconda3/envs/cdh2_r.env/lib/R/bin/Rscript\""
    exit 1
fi

filename="$1"
rscript_command="${2:-Rscript}"

if [ "$#" -gt 1 ]; then
    if [ ! -f "${rscript_command##* }" ]; then
        echo "ERROR: Nonexistent Rscript file path: ${rscript_command##* }"
        exit 1
    fi
fi

if [ ! -e "$filename" ]; then
    echo "ERROR: File not found: $filename"
    exit 1
fi

run_ipynb_cmd() {
	echo -e "\n\n\n---PROCESSING ${1}---\n\n\n"
	jupyter-nbconvert --to notebook --ExecutePreprocessor.timeout=-1 --execute $1
}

run_py_cmd() {
	echo -e "\n\n\n---PROCESSING ${1}---\n\n\n"
	python3.10 $1
}

run_sh_cmd() {
	echo -e "\n\n\n---PROCESSING ${1}---\n\n\n"
	bash $1
}

run_rmd_cmd() {
	echo -e "\n\n\n---RUNNING ${1}---\n\n\n"
	$rscript_command -e "rmarkdown::render('$1')"
}

run_env_cmd() {
    echo -e "\n\n\n---ENTERING ${1}---\n\n\n"
    source activate $1
}

run_denv_cmd() {
    echo -e "\n\n\n---EXITING CONDA ENV---\n\n\n"
    conda deactivate
}

run_cd_cmd() {
    dir_ext="$1"
    dir_no_ext="${dir_ext%.*}"
    echo -e "\n\n\n---CHANGING DIRECTORY TO ${dir_no_ext}---\n\n\n"
    cd $dir_no_ext
}

while IFS= read -r f_name; do
    # Skip empty lines and comments
    if [ -z "$f_name" ] || [[ "$f_name" == "#"* ]]; then
        continue
    fi

    # # Check for whitespace in the file name
    # if [[ "$f_name" =~ [[:space:]]+$ ]]; then
    #     echo "ERROR: File name contains whitespace: $f_name"
    #     exit 1
    # fi

    # Extract file extension
    extension="${f_name##*.}"
    # Remove anything after the extension
    extension="${extension%% *}"

    # # Check if the file exists and has the correct extension
    # if [ ! -e "$f_name" ] && [[ "$extension" != "env" ]] && [[ "$extension" != "denv" ]]; then
    #     echo "ERROR: File not found: $f_name"
    #     exit 1
    # fi

    # Determine the appropriate action based on the file extension
    case "$extension" in
        "ipynb")
            run_ipynb_cmd "$f_name"
            ;;
        "Rmd")
            run_rmd_cmd "$f_name"
            ;;
        "env")
            run_env_cmd "$f_name"
            ;;
        "denv")
            run_denv_cmd "$f_name"
            ;;
        "py")
            run_py_cmd "$f_name"
            ;;
        "sh")
            run_sh_cmd "$f_name"
            ;;
        "cd")
            run_cd_cmd "$f_name"
            ;;
        *)
            echo "Unknown extension: $extension for file $f_name"
            exit 1
            ;;
    esac

done < "$filename"

exit 0
