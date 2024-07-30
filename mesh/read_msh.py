import numpy as np


filename = "mesh/shape.msh"

try:
    with open(filename) as my_file:
        text = my_file.read()
    pass
except:
    print("error: unable to open mesh file!")
    exit(0)
print(f"reading {filename}..., done!")
