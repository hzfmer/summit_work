#!/usr/bin/zsh
set -e

base="../dhyp0.50_s1485839278_q100f00_orig_vs200"
list="dhyp0.50_s1485839278_q100f00_orig_vs500 dhyp0.50_s1485839278_q100f00_s05h005l100_vs200 dhyp1.00_s387100462_q100f00_orig_vs200 dhyp1.50_s372823598_q100f00_orig_vs200"

for dir in $list; do 
  break;
  echo $dir;
  cd $dir;
  ln -sf ${base}/IN3D.out .;
  ln -sf ${base}/param.sh .;
  ln -sf ${base}/stat.txt .;
  ln -sf ${base}/*slurm .;
  ln -sf /ccs/home/hzfmer/file_back/scripts/pickle_sites.py .;
  if [[ $dir != *"vs500"* ]]; then
    ln -sf ${base}/run.lsf .;
  fi;

  if [[ ! -f vel_sites.pickle ]]; then
    nstat=$(head -n 1 stat.txt)
    ndat=$(ls output_sfc/*dat | wc -l)
    if (( nstat != ndat )); then
      python pickle_sites.py ../stat_name_idx.txt 19440 14904 800;
    fi
  fi;

  # set to 0 > 0 if dont run slurms.
  if (( 0 > 0 )); then
    sbatch extrts.slurm;
    sbatch metrics.slurm;
    sbatch metrics2.slurm;
    sbatch metrics3.slurm;
    sbatch metrics4.slurm;
  fi
  cd ..;
done

for dir in $list; do
  cd ${dir};
  #rm extrts_31*;
  #sbatch extrts.slurm;
  python pickle_sites.py ../stat_name_idx.txt 19440 14904 800;
  cd ..;
done

exit -1;

list2="dhyp0.50_s1485839278_q100f00_s10h005l100_vs200 dhyp0.50_s1485839278_q100f00_s10h005l500_vs200 dhyp0.50_s1485839278_q100f00_s05h005l500_vs200"

for dir in $list2; do 
  echo $dir;
  cd $dir;
  ln -sf ${base}/IN3D.out .;
  ln -sf ${base}/param.sh .;
  ln -sf ${base}/stat.txt .;
  ln -sf ${base}/*slurm .;
  if [[ $dir != *"vs500"* ]]; then
    ln -sf ${base}/run.lsf .;
  fi;
  cd ..;
done
