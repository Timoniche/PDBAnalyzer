root_path="/Users/ddulaev/Documents/BioInformatics/SIMBA3D"
for (( i = 1; i < 22; i++))
do
  dir_path_iter="${root_path}/examples/command_line_tutorial/results/chr${i}"
  cd ${dir_path_iter}
  for file in *.pdb; do
    mv "$file" "Simba3d_chr${i}.pdb"
  done
done