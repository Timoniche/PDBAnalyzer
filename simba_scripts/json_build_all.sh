root_path="/Users/ddulaev/Documents/BioInformatics/SIMBA3D/examples/command_line_tutorial"
cd ${root_path}

for (( i = 4; i < 23; i++))
do
  dir_path="${root_path}/results/chr${i}"
  mkdir -p ${dir_path}
done

for (( i = 22; i < 23; i++))
do
  dir_path="${root_path}/results/chr${i}"
  mkdir -p dir_path
  json_task="tasks/taskChr${i}.json"
  simba3d -r ${json_task}
done
