import numpy as np

a = [1.348, 1.454, 1.381, 1.376, 1.453, 1.470]
b = 0.637
c = 0.652

for i in range(6):
    ans = np.sqrt(a[i]*a[i] - b*b + 0.5*c*c)
    print(f"layer{i+1} : {ans}")