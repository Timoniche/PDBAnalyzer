root_path="/Users/ddulaev/Documents/BioInformatics/SIMBA3D"
cd ${root_path}
for (( i = 5; i < 22; i++))
do
  dir_path_iter="${root_path}/examples/command_line_tutorial/results/chr${i}/*"
  for file in ${dir_path_iter}; do
    python simba3d_make_pdb.py ${file}
  done
done