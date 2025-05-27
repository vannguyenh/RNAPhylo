import os

model = input("Model: ")
MODEL_PATH = os.path.join('/Users/u7875558/Documents/RNAPhylo/allModels_SEED/outputs/', model)
# groups: raxml, raxmlP_wPseu, raxmlP_iPseu
GROUP = 'raxmlP_iPseu'

def control_output(path, group):
    control_path = os.path.join(path, group)
    notrun_rf = list()
    issue_rf = list()

    for rf in os.listdir(control_path):
        rf_path = os.path.join(control_path, rf)
        if os.path.isdir(rf_path):
            if len(os.listdir(rf_path)) == 10:
                print(f'RAxML cannot run {rf}')
                notrun_rf.append(rf)
            elif len(os.listdir(rf_path)) != 10 and len(os.listdir(rf_path)) != 50:
                print(f'{rf} has only {len(os.listdir(rf_path))} files.')
                issue_rf.append(rf)
            else:
                continue
        else:
            print(f'{rf} has issue with the pathway')
    return sorted(notrun_rf), issue_rf, len(os.listdir(control_path))

def main():
    notrunning_rf, problem_rf, len_rf = control_output(MODEL_PATH, GROUP)
    print(f'Not running RNAs: {notrunning_rf}')
    print(f'{problem_rf} is/are not done yet.')
    print(f'{len_rf} RNAs are inferred.')

if __name__ == "__main__":
    main()
