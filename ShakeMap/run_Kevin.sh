echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' records_pgv.txt)" > records_pgv.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 5 --intensity-file records_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/records;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' avg_pgv.txt)" > avg_pgv.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 5 --intensity-file avg_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/;
cd 223_5_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_5_18_uniform_pgv.txt)" > 223_5_18_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 5 --rupture-var-id 18 --intensity-file 223_5_18_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_5_18_uniform;
cd ..;
cd 223_6_22_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_22_uniform_pgv.txt)" > 223_6_22_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 6 --rupture-var-id 22 --intensity-file 223_6_22_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_22_uniform;
cd ..;
cd 223_11_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_26_uniform_pgv.txt)" > 223_11_26_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 11 --rupture-var-id 26 --intensity-file 223_11_26_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_26_uniform;
cd ..;
cd 223_10_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_25_uniform_pgv.txt)" > 223_10_25_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 10 --rupture-var-id 25 --intensity-file 223_10_25_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_25_uniform;
cd ..;
cd 223_10_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_26_uniform_pgv.txt)" > 223_10_26_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 10 --rupture-var-id 26 --intensity-file 223_10_26_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_26_uniform;
cd ..;
cd 223_4_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_18_uniform_pgv.txt)" > 223_4_18_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 4 --rupture-var-id 18 --intensity-file 223_4_18_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_18_uniform;
cd ..;
cd 223_11_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_25_uniform_pgv.txt)" > 223_11_25_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 11 --rupture-var-id 25 --intensity-file 223_11_25_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_25_uniform;
cd ..;
cd 223_6_17_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_17_uniform_pgv.txt)" > 223_6_17_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 6 --rupture-var-id 17 --intensity-file 223_6_17_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_17_uniform;
cd ..;
cd 223_4_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_19_uniform_pgv.txt)" > 223_4_19_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 4 --rupture-var-id 19 --intensity-file 223_4_19_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_19_uniform;
cd ..;
cd 223_3_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_3_19_uniform_pgv.txt)" > 223_3_19_uniform_pgv.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period -1 --source-id 223 --rupture-id 3 --rupture-var-id 19 --intensity-file 223_3_19_uniform_pgv.txt --spacing 0.005 --colorbar-min 0.5386623740196228 --colorbar-max 175.64414978027344 -o /Users/zhh076/work/ShakeMap/CyberShake/223_3_19_uniform;
cd ..;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' records_pga.txt)" > records_pga.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 5 --intensity-file records_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/records;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' avg_pga.txt)" > avg_pga.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 5 --intensity-file avg_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/;
cd 223_5_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_5_18_uniform_pga.txt)" > 223_5_18_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 5 --rupture-var-id 18 --intensity-file 223_5_18_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_5_18_uniform;
cd ..;
cd 223_6_22_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_22_uniform_pga.txt)" > 223_6_22_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 6 --rupture-var-id 22 --intensity-file 223_6_22_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_22_uniform;
cd ..;
cd 223_11_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_26_uniform_pga.txt)" > 223_11_26_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 11 --rupture-var-id 26 --intensity-file 223_11_26_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_26_uniform;
cd ..;
cd 223_10_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_25_uniform_pga.txt)" > 223_10_25_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 10 --rupture-var-id 25 --intensity-file 223_10_25_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_25_uniform;
cd ..;
cd 223_10_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_26_uniform_pga.txt)" > 223_10_26_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 10 --rupture-var-id 26 --intensity-file 223_10_26_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_26_uniform;
cd ..;
cd 223_4_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_18_uniform_pga.txt)" > 223_4_18_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 4 --rupture-var-id 18 --intensity-file 223_4_18_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_18_uniform;
cd ..;
cd 223_11_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_25_uniform_pga.txt)" > 223_11_25_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 11 --rupture-var-id 25 --intensity-file 223_11_25_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_25_uniform;
cd ..;
cd 223_6_17_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_17_uniform_pga.txt)" > 223_6_17_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 6 --rupture-var-id 17 --intensity-file 223_6_17_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_17_uniform;
cd ..;
cd 223_4_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_19_uniform_pga.txt)" > 223_4_19_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 4 --rupture-var-id 19 --intensity-file 223_4_19_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_19_uniform;
cd ..;
cd 223_3_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_3_19_uniform_pga.txt)" > 223_3_19_uniform_pga.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0 --source-id 223 --rupture-id 3 --rupture-var-id 19 --intensity-file 223_3_19_uniform_pga.txt --spacing 0.005 --colorbar-min 0.00726829981431365 --colorbar-max 1.1820369958877563 -o /Users/zhh076/work/ShakeMap/CyberShake/223_3_19_uniform;
cd ..;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' avg_psa05.txt)" > avg_psa05.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 5 --intensity-file avg_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/;
cd 223_5_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_5_18_uniform_psa05.txt)" > 223_5_18_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 5 --rupture-var-id 18 --intensity-file 223_5_18_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_5_18_uniform;
cd ..;
cd 223_6_22_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_22_uniform_psa05.txt)" > 223_6_22_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 6 --rupture-var-id 22 --intensity-file 223_6_22_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_22_uniform;
cd ..;
cd 223_11_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_26_uniform_psa05.txt)" > 223_11_26_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 11 --rupture-var-id 26 --intensity-file 223_11_26_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_26_uniform;
cd ..;
cd 223_10_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_25_uniform_psa05.txt)" > 223_10_25_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 10 --rupture-var-id 25 --intensity-file 223_10_25_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_25_uniform;
cd ..;
cd 223_10_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_26_uniform_psa05.txt)" > 223_10_26_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 10 --rupture-var-id 26 --intensity-file 223_10_26_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_26_uniform;
cd ..;
cd 223_4_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_18_uniform_psa05.txt)" > 223_4_18_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 4 --rupture-var-id 18 --intensity-file 223_4_18_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_18_uniform;
cd ..;
cd 223_11_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_25_uniform_psa05.txt)" > 223_11_25_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 11 --rupture-var-id 25 --intensity-file 223_11_25_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_25_uniform;
cd ..;
cd 223_6_17_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_17_uniform_psa05.txt)" > 223_6_17_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 6 --rupture-var-id 17 --intensity-file 223_6_17_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_17_uniform;
cd ..;
cd 223_4_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_19_uniform_psa05.txt)" > 223_4_19_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 4 --rupture-var-id 19 --intensity-file 223_4_19_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_19_uniform;
cd ..;
cd 223_3_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_3_19_uniform_psa05.txt)" > 223_3_19_uniform_psa05.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 0.5 --source-id 223 --rupture-id 3 --rupture-var-id 19 --intensity-file 223_3_19_uniform_psa05.txt --spacing 0.005 --colorbar-min 0.0137948 --colorbar-max 2.0595357 -o /Users/zhh076/work/ShakeMap/CyberShake/223_3_19_uniform;
cd ..;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' records_psa10.txt)" > records_psa10.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 5 --intensity-file records_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/records;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' avg_psa10.txt)" > avg_psa10.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 5 --intensity-file avg_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/;
cd 223_5_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_5_18_uniform_psa10.txt)" > 223_5_18_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 5 --rupture-var-id 18 --intensity-file 223_5_18_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_5_18_uniform;
cd ..;
cd 223_6_22_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_22_uniform_psa10.txt)" > 223_6_22_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 6 --rupture-var-id 22 --intensity-file 223_6_22_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_22_uniform;
cd ..;
cd 223_11_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_26_uniform_psa10.txt)" > 223_11_26_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 11 --rupture-var-id 26 --intensity-file 223_11_26_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_26_uniform;
cd ..;
cd 223_10_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_25_uniform_psa10.txt)" > 223_10_25_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 10 --rupture-var-id 25 --intensity-file 223_10_25_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_25_uniform;
cd ..;
cd 223_10_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_26_uniform_psa10.txt)" > 223_10_26_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 10 --rupture-var-id 26 --intensity-file 223_10_26_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_26_uniform;
cd ..;
cd 223_4_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_18_uniform_psa10.txt)" > 223_4_18_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 4 --rupture-var-id 18 --intensity-file 223_4_18_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_18_uniform;
cd ..;
cd 223_11_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_25_uniform_psa10.txt)" > 223_11_25_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 11 --rupture-var-id 25 --intensity-file 223_11_25_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_25_uniform;
cd ..;
cd 223_6_17_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_17_uniform_psa10.txt)" > 223_6_17_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 6 --rupture-var-id 17 --intensity-file 223_6_17_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_17_uniform;
cd ..;
cd 223_4_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_19_uniform_psa10.txt)" > 223_4_19_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 4 --rupture-var-id 19 --intensity-file 223_4_19_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_19_uniform;
cd ..;
cd 223_3_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_3_19_uniform_psa10.txt)" > 223_3_19_uniform_psa10.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 1 --source-id 223 --rupture-id 3 --rupture-var-id 19 --intensity-file 223_3_19_uniform_psa10.txt --spacing 0.005 --colorbar-min 0.0053109 --colorbar-max 1.2342584 -o /Users/zhh076/work/ShakeMap/CyberShake/223_3_19_uniform;
cd ..;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' records_psa30.txt)" > records_psa30.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 5 --intensity-file records_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/records;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' avg_psa30.txt)" > avg_psa30.txt;
java -Xmx2G -cp opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 5 --intensity-file avg_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/;
cd 223_5_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_5_18_uniform_psa30.txt)" > 223_5_18_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 5 --rupture-var-id 18 --intensity-file 223_5_18_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_5_18_uniform;
cd ..;
cd 223_6_22_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_22_uniform_psa30.txt)" > 223_6_22_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 6 --rupture-var-id 22 --intensity-file 223_6_22_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_22_uniform;
cd ..;
cd 223_11_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_26_uniform_psa30.txt)" > 223_11_26_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 11 --rupture-var-id 26 --intensity-file 223_11_26_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_26_uniform;
cd ..;
cd 223_10_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_25_uniform_psa30.txt)" > 223_10_25_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 10 --rupture-var-id 25 --intensity-file 223_10_25_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_25_uniform;
cd ..;
cd 223_10_26_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_10_26_uniform_psa30.txt)" > 223_10_26_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 10 --rupture-var-id 26 --intensity-file 223_10_26_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_10_26_uniform;
cd ..;
cd 223_4_18_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_18_uniform_psa30.txt)" > 223_4_18_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 4 --rupture-var-id 18 --intensity-file 223_4_18_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_18_uniform;
cd ..;
cd 223_11_25_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_11_25_uniform_psa30.txt)" > 223_11_25_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 11 --rupture-var-id 25 --intensity-file 223_11_25_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_11_25_uniform;
cd ..;
cd 223_6_17_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_6_17_uniform_psa30.txt)" > 223_6_17_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 6 --rupture-var-id 17 --intensity-file 223_6_17_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_6_17_uniform;
cd ..;
cd 223_4_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_4_19_uniform_psa30.txt)" > 223_4_19_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 4 --rupture-var-id 19 --intensity-file 223_4_19_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_4_19_uniform;
cd ..;
cd 223_3_19_uniform;
echo "$(awk '{if($1 < 0) {t = $1; $1 = $2; $2=t}; print}' 223_3_19_uniform_psa30.txt)" > 223_3_19_uniform_psa30.txt;
java -Xmx2G -cp ../opensha-cybershake-all.jar org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator --study STUDY_15_12 --period 3 --source-id 223 --rupture-id 3 --rupture-var-id 19 --intensity-file 223_3_19_uniform_psa30.txt --spacing 0.005 --colorbar-min 0.0005351 --colorbar-max 0.2636306 -o /Users/zhh076/work/ShakeMap/CyberShake/223_3_19_uniform;
cd ..;
