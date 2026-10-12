import os
import parameters
import hashing
from USER_VARIABLES import RESULTS_DIRECTORY

ids = set()
cpids = set()
for result_dir in os.listdir(RESULTS_DIRECTORY):
    rd = os.path.join(RESULTS_DIRECTORY, result_dir)
    if not os.path.isdir(rd): continue
        
    for f in os.listdir(rd):
        if f.startswith('.'): continue
        fp = os.path.join(rd, f)
        try:
            p, l, r = parameters.load_parameter_file(fp)
            model_id = hashing.build_model_id(r['model'])
            cp_mid = [_k for _k in f.split('_')
                        if _k.startswith('model-')][0]
            cpids.add(cp_mid)
            print(model_id, f)
            ids.add(model_id)
            break
        except Exception as ex:
            print(ex)
            input()
            continue
print(ids)
print(cpids)
