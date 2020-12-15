import numpy as np
import subprocess



exe_file = "../main.exe"

T = 0.5
DT = [5e-1, 5e-2, 5e-3, 5e-4, 5e-5, 5e-6]

for dt in DT:
    print(dt)
    subprocess.call([exe_file, str(dt), str(T)])
