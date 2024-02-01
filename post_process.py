from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import numpy as np


def post_process(x):
    [rel_depth, t1, t2, t3] = x
    project = "waterbomb"
    report_address = "Report.txt"

    odb_name = project + '.odb'
    odb = session.openOdb(name=odb_name)
    data = session.XYDataFromHistory(name='data', odb=odb, steps=('Analysis',),
                                     outputVariableName='Strain energy: ALLSE for Whole Model')
    energy = np.array(data)[:, 1]
    displacement = np.array([odb.steps['Analysis'].frames[i].frameValue for i in range(len(energy))])
    derivative = np.array(
        [(energy[i + 1] - energy[i]) / (displacement[i + 1] - displacement[i]) for i in range(len(energy) - 2)])

    end = derivative.tolist().index(derivative[derivative < 0][-1])
    start = end - 1
    iteration = derivative[start]
    while iteration < 0:
        start -= 1
        iteration = derivative[start]

    if (len(derivative[derivative < 0]) > 0) and (float(start) / len(derivative) > 0.25):
        delta = (energy[start + 2] - energy[end + 1]) / energy[start + 2]
    else:
        delta = -float('inf')

    with open(report_address, 'a') as fp:
        fp.write(str(rel_depth) + "\t" +
                 str(t1) + "\t" +
                 str(t2) + "\t" +
                 str(t3) + "\t" +
                 str(delta) + "\n")

    odb.close()
    return
