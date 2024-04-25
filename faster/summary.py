import os
import sys
import USER_VARIABLES
import model
import data
from datetime import datetime

def summary(args):
    d = data.get_data_before_day()
    results = {}
    result_source = USER_VARIABLES.LOG_DIRECTORY
    for f in os.listdir(result_source):
        parameter_source = os.path.join(result_source, f)
        if not os.path.isdir(parameter_source) or not f.startswith('fit'):
            continue
        if len(os.listdir(parameter_source)) == 0:
            continue
        if args and not all([a in f for a in args]):
            continue
        best_loss, best_parameters = model.get_best_loss_parameters(parameter_source)
        model_type = 'complex' if 'Fe3_v_max' in best_parameters else 'simple'
        loaded_model = model.Model(model.get_pathways(model_type))
        loaded_model.parameters().set(best_parameters)
 
        clean_name = f.replace('fit_','')
        clean_name = clean_name.replace('simple','').replace('complex','')
        clean_name = clean_name.replace('log', '')
        clean_name, date = clean_name.split('2024')
        replicas = [s for s in clean_name.split('_') if not s == '']
        runs = {replica_name: None for replica_name in replicas}

        for replica_name in replicas:
            replica = d[replica_name]
            run = loaded_model.predict(replica)
            runs[replica_name] = run
        
        results[f] = {'loss': best_loss,
                      'parameters': best_parameters,
                      'date': date,
                      'R2': {r_name: r['R2'] for r_name, r in runs.items()}}
        
    return results

if __name__ == '__main__':
    args = sys.argv[1:]
    only_today = False
    if 'today' in args:
        only_today = True
        args.remove('today')
    res = summary(args)
    for k in sorted(res.keys()):
        v = res[k]
        best_loss = v['loss']
        
        date = datetime.strptime('2024' + v['date'], '%Y-%m-%d--%H-%M-%S')
        if only_today and date < datetime.today().replace(hour = 0, minute = 0, second = 0):
            continue
        st = f'{k[:55]:<55s}  {best_loss:8.4f}'
        print()
        print(st)
        for r_name, r2 in v['R2'].items():
            r2_str = f'   R2 (CO2) = {r2["CO2"]:5.2f} ({r_name})\n   R2 (CH4) = {r2["CH4"]:5.2f} ({r_name})'
            print(r2_str)
 
