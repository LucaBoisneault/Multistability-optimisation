import PyNomad
from functions import *

name = "waterbomb"

x0 = [1.0, 0.1, 0.5, 0.9]

lb = [0.0, 0.01, 0.01, 0.01]
ub = [2.0, 0.99, 0.99, 0.99]

params = ["DIMENSION "+str(int(len(x0))),
          "BB_OUTPUT_TYPE OBJ", 
          "MAX_BB_EVAL 1000", 
          "BB_INPUT_TYPE (R R R R)",
          "DISPLAY_ALL_EVAL true",
          "DISPLAY_STATS BBE ( sol ) OBJ", 
          'stats_file "Blackbox_result.txt" BBE ( sol ) OBJ',
          'DISPLAY_DEGREE 2']

PyNomad.optimize(bb_pynomad, x0, lb, ub, params)
