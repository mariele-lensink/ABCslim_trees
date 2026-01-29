#command line scripts for running each model

#bottleneck
slim \
  -d "ID='TESTBOT'" \
  -d "OUTTREES='trees/bottleneck'" \
  -d END=200 \
  -d N_END=1135 \
  -d SELF=0.98 \
  -d T_BOT=75 -d T_REC=175 -d BOT_FRAC=0.2 \
  -d gmu=1e-8 -d imu=1e-8 -d gd=0.3 -d id=0.3 -d gdfe=-0.01 -d idfe=-0.01 \
  ABCtrees_bottleneck.slim

# Convert bottleneck trees to VCF
python3 trees_to_vcf.py \
  --trees trees/bottleneck/TESTBOT.trees \
  --out_vcf trees/bottleneck/TESTBOT.vcf \
  --gene_csv genomeinfo/gene100.csv \
  --intergene_csv genomeinfo/intergene100.csv \
  --L 1000000 --gmu 1e-8 --imu 1e-8 --gd 0.3 --id 0.3

#expansion
slim \
  -d "ID='TESTEXP'" \
  -d "OUTTREES='trees/expansion'" \
  -d END=200 \
  -d N_END=1135 \
  -d SELF=0.98 \
  -d T_EXP=150 -d EXP_MULT=2.0 \
  -d gmu=1e-8 -d imu=1e-8 -d gd=0.3 -d id=0.3 -d gdfe=-0.01 -d idfe=-0.01 \
  ABCtrees_expansion.slim

# Convert expansion trees to VCF
python3 trees_to_vcf.py \
  --trees trees/expansion/TESTEXP.trees \
  --out_vcf trees/expansion/TESTEXP.vcf \
  --gene_csv genomeinfo/gene100.csv \
  --intergene_csv genomeinfo/intergene100.csv \
  --L 1000000 --gmu 1e-8 --imu 1e-8 --gd 0.3 --id 0.3

#constant size
slim \
 -d "ID='TEST001'" \
 -d "OUTTREES='trees/constant'" \
 -d END=200   -d SIMPLIFY_INTERVAL=50 \
 -d NPOP=50   \
 -d SELF=0.98   \
 -d gmu=1e-8   -d imu=1e-8   -d gd=0.3   -d id=0.3   -d gdfe=-0.01   -d idfe=-0.01\
 ABCtrees.slim 2>&1 | tee trees/constant/slim_test_run.log

# Convert constant size trees to VCF
python3 trees_to_vcf.py \
  --trees trees/constant/TEST001.trees \
  --out_vcf trees/constant/TEST001.vcf \
  --gene_csv genomeinfo/gene100.csv \
  --intergene_csv genomeinfo/intergene100.csv \
  --L 1000000 --gmu 1e-8 --imu 1e-8 --gd 0.3 --id 0.3