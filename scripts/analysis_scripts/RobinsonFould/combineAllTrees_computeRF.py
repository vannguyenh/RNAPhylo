import os
import subprocess

def produceCombinedTrees(dir_input, dir_output):
    os.makedirs(dir_output, exist_ok=True)
    # endings for output files
    endings={
        "raxml": "raxml",
        "raxmlP_wPseu": "raxmlPw",
        "raxmlP_iPseu" : "raxmlPi"
    }

    # List all files in the input directory
    methods = list(endings.keys())
    rnas = os.listdir(os.path.join(dir_input, methods[0]))
    for method in methods:
        for rna in rnas:
            input_method = os.path.join(dir_input, method, rna)
            input_files = sorted([f for f in os.listdir(input_method) if f.startswith("RAxML_bestTree.")])
            output_rna=os.path.join(dir_output, rna)
            os.makedirs(output_rna, exist_ok=True)
            output_file=os.path.join(output_rna, f"{rna}.{endings[method]}")
            with open(output_file, "w") as outfile:
                for filename in input_files:
                    file_path = os.path.join(input_method, filename)
                    with open(file_path, "r") as infile:
                        content = infile.read()
                        outfile.write(content)

def run_command(command):
    process=subprocess.Popen(command, shell=True)
    process.communicate()

def computeRFdistance_iqtreecmd(dcombine_path):
    rnas=os.listdir(dcombine_path)
    for rna in rnas:
        dir_combine_rna=os.path.join(dcombine_path, rna)
        for f in os.listdir(dir_combine_rna):
            if f.endswith("raxml"):
                raxTree=os.path.join(dir_combine_rna, f)
            elif f.endswith("raxmlPw"):
                raxPwPTree=os.path.join(dir_combine_rna, f)
            else:
                raxPwoPTree=os.path.join(dir_combine_rna, f)
        #prefix=f"{dcombine_path}/{rna}/"
        command=f"bash computeRFdistance.sh {raxTree} {raxPwPTree} {raxPwoPTree} {dir_combine_rna} {rna}"
        run_command(command)

def main():
    # Directory containing the input files
    DIR_OUTPUTS = "//Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/outputs"
    MODEL = input('Model: ')
    dir_combined = os.path.join(DIR_OUTPUTS, 'Robinson_Foulds', MODEL)

    produceCombinedTrees(os.path.join(DIR_OUTPUTS, MODEL), dir_combined)
    computeRFdistance_iqtreecmd(dir_combined)

if __name__=="__main__":
    main()


