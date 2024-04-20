import os

import USER_VARIABLES
import model
import data

def summary():
    d = data.get_data_before_day()
    results = {}
    result_source = USER_VARIABLES.LOG_DIRECTORY
    for f in os.listdir(result_source):
        parameter_source = os.path.join(result_source, f)
        if not os.path.isdir(parameter_source) or not f.startswith('fit'):
            continue

        best_loss, best_parameters = model.get_best_loss_parameters(parameter_source)
        model_type = 'complex' if 'M_Fe3' in best_parameters else 'simple'
        loaded_model = model.Model(model.get_pathways(model_type))
        loaded_model.parameters().set(best_parameters)
  
        replicas = [s for s in f.replace('fit_','').replace('log', '').split('2024')[0].split('_') if not s == '']
        runs = {replica_name: None for replica_name in replicas}
        for replica_name in replicas:
            replica = d[replica_name]
            run = loaded_model.predict(replica)
            runs[replica_name] = run
        
        results[f] = {'loss': best_loss,
                      'parameters': best_parameters,
                      'R2': {r_name: r['R2'] for r_name,r in runs.items()}}
        
    return results

if __name__ == '__main__':
    res = summary()
    for k, v in res.items():
        print()
        best_loss = v['loss']
        st = f'{k[:25]:<25s}  {best_loss:8.4f}'
        print(st)
        for r_name, r2 in v['R2'].items():
            r2_str = f'   R2 (CO2) = {r2["CO2"]:5.2f} ({r_name})\n   R2 (CH4) = {r2["CH4"]:5.2f} ({r_name})'
            print(r2_str)
 
