#! /bin/bash
#SBATCH --job-name=zfc
#SBATCH --output=/gpfs/share/home/1801111526/job.out/%x.%j_%A_%a.%N.out
#SBATCH --error=/gpfs/share/home/1801111526/job.out/%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0256G
#SBATCH --qos=normal
#SBATCH --get-user-env
#SBATCH -n 20 
#SBATCH --cpu-freq=high
#SBATCH -A hpc0006169510
#SBATCH --mail-type=end
#SBATCH --mail-user=1801111526@pku.edu.cn
#SBATCH --time=120:00:00

#path='/gpfs/share/home/1801111526/project/sty_ptm/downsample'
#for i in lib1_0.2 lib1_0.3 lib1_0.5 lib1_0.8 lib2_0.2 lib2_0.3 lib2_0.5 lib2_0.8
#do
#./zfc -i ${path}/${i}_input.txt -o ${path}/${i} --plot
#done

./zfc -i /gpfs/share/home/1801111526/project/crc/screen/cbe/drug2/dmso2_d0_bbk_zfc.txt -o /gpfs/share/home/1801111526/project/crc/screen/cbe/drug2/drug2_bbk_zfc/dmso2_bbk --plot
