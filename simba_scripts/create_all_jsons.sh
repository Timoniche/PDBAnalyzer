for (( i = 1; i < 23; i++))
do
  file_path="/Users/ddulaev/Documents/BioInformatics/SIMBA3D/examples/command_line_tutorial/tasks/taskChr${i}.json"
  body="[[{
    \"file_names_inputdir\":\"data/\",
    \"file_names_pairwise_contact_matrix\":\"chr${i}.npy\",
    \"file_names_outputdir\":\"results/chr${i}\"
    }]]"
  echo ${body} > ${file_path}
done