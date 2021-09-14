#script to run the Pipeline multiple times with the specified config files in the "/config_files" directory
for f in trnl_CH ITS2_config
do
  bash run_pipeline.sh $PWD/config_files/$f.yaml
done
