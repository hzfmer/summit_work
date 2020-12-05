#!/usr/bin/zsh

StringVal="q50f00_orig q100f00_orig q100f04_orig q100f06_orig q100f08_orig q100f00_s05h005l100 q100f04_s05h005l100 q100f06_s05h005l100 q100f06_s05h005l200 q100f06_s05h005l500 q100f06_s05h010l100 q100f06_s05h010l500 q100f06_s05h015l100 q100f06_s10h005l100 q100f08_s05h005l100 q50f00_s05h005l100 q100f00_orig_vs500"
StringVal="q50f00_orig q100f00_orig q50f00_s05h005l100 q100f00_s05h005l100 q100f00_orig_vs500 q100f04_orig q100f06_orig q100f08_orig q100f04_s05h005l100 q100f06_s05h005l100 q100f06_s05h005l200 q100f06_s05h005l500 q100f06_s05h010l100 q100f06_s05h010l500 q100f06_s05h015l100 q100f06_s10h005l100 q100f08_s05h005l100"


for val in $StringVal; do
#  continue
  echo $val;
  cd $val;
  echo $PWD;
#  ls -l mesh_0
#  ls -l mesh_1
#  ls -l source_0
#  ls -l source_1
#  bsub run_origcode.lsf
#  tail -5 param_origcode.sh
#  find ./ -name "*lock" -exec rm {} +;
#  rm 2*{err,out}
#  ln -sf ../stat_small.txt stat.txt;
#  ln -sf ~/file_back/scripts/pickle_sites.py .;
#  ln -sf ~/file_back/scripts/generate_in3dout.py .;
#  if [[ "" && ($val == *"topo"* || $val == *"vs500"*) ]]; then
#    cp ../topo_qf06_s05h005l100/IN3D.out .;
#  else
#    cp ../noqf_orig/IN3D.out .;
#  fi
#  find ./ -mtime +1 -name "*bin" -exec rename Hz_ _ {} +
#   find ./ -name "*bin" -exec rename Hz_ _ {} +
#  cp ../topo_qf06_s05h005l100/metrics.slurm .;
#  sbatch metrics.slurm;
  #  sbatch gmrot.slurm;
#  sleep 10m;
#  if [[ "" && ! -f "gmrotD50_00.50Hz.bin" ]]; then 
#    echo "yes";
#    cp ../topo_qf06_s05h005l100/gmrot.slurm .;
#    sleep 8m;
#  fi
#  cp ../noqf_orig/IN3D.out .;
#  cp ../topo_qf06_s05h005l100/extrts.slurm .;
#  cp ../topo_qf06_s05h005l100/gmrot.slurm .;
   diff ../q100f00_orig/metrics.slurm ./metrics.slurm;
#  cp ../topo_qf06_s05h005l100/temp_metrics.slurm .;
#  sbatch metrics.slurm;
#  if [[ "" && -f vel_sites.pickle ]]; then
#    rm output_sfc/S*dat;
#    rm output_sfc/S*_1_*;
#  fi
  if [[ -n $(find ./ -name "pgv*bin" -size -47775744c) ]] ; then 
    echo yes; 
#    sbatch temp_metrics.slurm;
  fi
#  python pickle_sites.py ../la_habra_small_statlist_3456.idx 3456 3456;
  if [[ $val == "noqf_orig" || $val == "topo_noqf_orig" ]]; then
    :
  fi
  if [[ $val == *"topo"* ]]; then
    #bsub run.lsf;
    echo "Check topo"
  else
    #bsub run_origcode.lsf;
    echo "Not topo"
  fi
#  if [[ ! -n $(find ./ -name "pgv*bin" -mtime -1) ]]; then
#    sbatch metrics.slurm;
#  fi
#  if [[ ! -n $(find ./ -name "gmrot*bin" -mtime -1) ]]; then
#    sbatch gmrot.slurm;
#  fi
#  if [[ ! -n $(find ./output_sfc/ -name "*dat" -mtime 1) ]]; then
#    sbatch extrts.slurm;
#  fi
  # sleep 15m;
  cd ..;
done

#StringVal="topo_q50f00_s05h005l050 topo_q50f00_s05h005l100 topo_q50f00_s05h005l200 topo_q50f00_s05h005l500 topo_q50f00_s05h010l100 topo_q50f00_s05h015l100 topo_q50f00_s10h005l100 topo_q50f00_s75h005l100 topo_q50f02_s05h005l100 topo_q50f04_s05h005l100 topo_q50f06_s05h005l100 topo_q50f08_s05h005l100 topo_q100f00_s05h005l100 topo_q100f02_s05h005l100 topo_q100f04_s05h005l100 topo_q100f08_s05h005l100"
StringVal="topo_q50f02_orig topo_q50f04_orig topo_q50f06_orig topo_q50f08_orig topo_q100f02_orig topo_q100f04_orig topo_q100f06_orig topo_q100f08_orig"
for val in $StringVal; do
  continue
  cd $val;
  echo $PWD;
#  sbatch extrts.slurm;
#  sbatch gmrot.slurm;
#  sbatch metrics.slurm;
#  python pickle_sites.py 3456 3456 ../la_habra_small_statlist_3456.idx;
  if [[ -n $(find -name "S._1_*") ]]; then
    echo "Exist S_1_* file"
    rm output_sfc/S*_1_*000
  fi
  cd ..;
done

exit 0;
