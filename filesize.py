import os
import pathlib
path = pathlib.Path().resolve()

filesizesComplete = []
filesizesRay = []

for filepath in os.listdir():
    if filepath[-2] == 'h':
        filesizesRay.append(os.path.getsize(filepath))
    elif filepath[-2] == 'f':
        filesizesComplete.append(os.path.getsize(filepath))

import numpy as np
import matplotlib.pyplot as plt

print(np.mean(filesizesComplete))
print(np.mean(filesizesRay))
print(np.std(filesizesComplete))
print(np.std(filesizesRay))
