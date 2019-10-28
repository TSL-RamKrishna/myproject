for script in  scripts/*transcripts.sh; do sbatch --mem 20G -o log/$(basename  $script | sed 's/.sh$/.log/') -J $(basename $script | sed 's/.sh$//') $script;  done
