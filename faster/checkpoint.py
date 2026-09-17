from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

# RUN hash
# model hash
# 

# configure callback:
#   keep only N best, if N is None, keep all
#   store location

class PrintCallback():
    def __call__(self, objective):
        args, kwargs, loss_value = objective.last_call()
        _, _, best_loss = objective.best_call()
        cnt = objective.call_count()
        print(f'\rcall {cnt:6d}: loss value {loss_value:8.3g}, best loss: {best_loss:8.3g}', end = '')

class CheckpointCallback():
    def __call__(self, objective):
        pass