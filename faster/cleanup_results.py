# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 09:45:31 2026

@author: Lara
"""
import os
import USER_VARIABLES

target = USER_VARIABLES.LOG_DIRECTORY

all_folders = []

for d in os.listdir(target):
    folder = os.path.join(target, d)
    
    all_folders.append((len(os.listdir(folder)), folder))
                       

all_folders = sorted(all_folders, reverse=True)

for _, f in all_folders:
    print(f)
    files = os.listdir(f)
    calls = [int(f.split('_')[1]) for f in files]
    
    p = sorted(zip(calls, files))
    q = list(zip(*p))[1]
    best_file = os.path.join(f, q[-1])
    
    y = input(str(len(files)) + ' files, best: ' + q[-1])

    if y == 'y':
        print('deleting')
        for file in files:
            path = os.path.join(f,file)
            if path == best_file: 
                continue
            os.remove(path)
