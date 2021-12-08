#!/bin/sh

# Execute the job in current working directory. All inputs and outputs will be directed from/to your present working directory.
#$ -cwd

# Keep my current environment variables
#$ -V

# Set the name of the job
#$ -N example

# Set output and error files.
#$ -o ./results/coso.out
#$ -e ./results/coso.error

# Choose the parallel environment and set how many cores you need. (Only necessary for parallel jobs.)
##$ -pe MPI 8

# Choose your queue and host/hostgroup. For more info please see https://www.engr.colostate.edu/ens/info/researchcomputing/cluster/queueinfo.html
##$ -q defaultfaculty.q@@students
##$ -q short.q
#$ -q jathar.q

echo STARTED @ $(date)

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/students/yicongh/Programs/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/students/yicongh/Programs/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/students/yicongh/Programs/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/students/yicongh/Programs/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

python runme.py 0 SOM-TOMAS

echo ENDED @ $(date)

