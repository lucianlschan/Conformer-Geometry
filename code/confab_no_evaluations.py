import re, os
import numpy as np
import pandas as pd

# Get the number of conformers evaluated in Confab

def check_log(file_in):
    name = []
    target= []
    for paths, dirs, files in os.walk(file_in):
        if paths == file_in:
            continue
        for x in files:
            if x.endswith('log.txt'):
                f = open(os.path.join(paths, x), 'r')
                for l in f.readlines():
                    match = re.match("..generated (\d+) conformers", l)
                    if match: name.append(int (match.groups()[0]))
                f.close()
        target.append(os.path.basename(paths))
    return name, target

source_1 = "/path/to/data/"
result = check_log(source_1)
data_frame = pd.DataFrame({"target": result[1], "Number": result[0]})
print(np.max(result[0]))
print(result[1][np.argmax(result[0])])
data_frame.to_csv(os.path.dirname(os.getcwd()) + "/result/confab_no_evaluation/six_confab_no.txt", index=False)
