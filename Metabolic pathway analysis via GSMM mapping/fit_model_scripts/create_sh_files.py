for i in range(150):
        with open(f"submit2cluster_{i}.sh", "w") as file:
                file.write("#!/bin/bash")
                file.write("\n")
                file.write(f"Rscript --vanilla fit_reaction_models.R -i abund_and_metadata.csv -c rxn_names.csv -r FALSE -b {i*51+2} -e {min((i+1)*51+1, 7605)}")