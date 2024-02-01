from abaqus import *
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


def model(x):
    [rel_depth, t1, t2, t3] = x
    project = 'waterbomb'

    model_name = project + '_model'
    cae_file = project + '.cae'
    odb_name = project + '.odb'

    # --------------------------------------------------------------------
    # Fixed variable that could be implemented in the optimization, but would result in a divergence
    # --------------------------------------------------------------------

    n = 4
    r_ext = 10.0
    r_int = 1.0
    size = 0.25
    size2 = 0.1

    quadratic = False
    hyperelastic = False

    hard_mat = 2410.0
    soft_mat = 15.2
    c_coefficient = 0.015

    m = mdb.Model(modelType=STANDARD_EXPLICIT, name=model_name)

    # --------------------------------------------------------------------
    # part 1 - base faces
    # --------------------------------------------------------------------

    s = m.ConstrainedSketch(name='base', sheetSize=200.0)
    s.ArcByCenterEnds(center=(0.0, 0.0), direction=COUNTERCLOCKWISE,
                      point1=(r_int, 0.0),
                      point2=(r_int * np.cos(np.pi / n), r_int * np.sin(np.pi / n)))

    s.ArcByCenterEnds(center=(0.0, 0.0), direction=COUNTERCLOCKWISE,
                      point1=(r_ext, 0.0),
                      point2=(r_ext * np.cos(np.pi / n), r_ext * np.sin(np.pi / n)))

    s.Line(point1=(r_int, 0.0),
           point2=(r_ext, 0.0))

    s.Line(point1=(r_int * np.cos(np.pi / n), r_int * np.sin(np.pi / n)),
           point2=(r_ext * np.cos(np.pi / n), r_ext * np.sin(np.pi / n)))

    p = m.Part(dimensionality=THREE_D, name='base', type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s)

    p.DatumAxisByPrincipalAxis(principalAxis=YAXIS)

    s1 = m.ConstrainedSketch(gridSpacing=1.0, name='__profile__', sheetSize=100.0, transform=p.MakeSketchTransform(
        sketchPlane=p.faces[0],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=p.datums[2],
        sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))

    s1.Line(point1=(r_ext * np.cos(t1 * np.pi / n), r_ext * np.sin(t1 * np.pi / n)),
            point2=(r_int * np.cos(t2 * np.pi / n), r_int * np.sin(t2 * np.pi / n)))

    s1.Line(point1=(r_int * np.cos(t2 * np.pi / n), r_int * np.sin(t2 * np.pi / n)),
            point2=(r_ext * np.cos(t3 * np.pi / n), r_ext * np.sin(t3 * np.pi / n)))

    s1.ConstructionLine(point1=(0.0, 0.0), point2=(np.cos(np.pi / n), np.sin(np.pi / n)))
    p.PartitionFaceBySketch(faces=p.faces[0], sketch=s1,
                            sketchUpEdge=p.datums[2])

    p.Set(faces=(p.faces[0:2]), name='soft')
    p.Set(faces=(p.faces[:]), name='all')
    p.Set(faces=(p.faces[2:3]), name='hard')
    p.Set(vertices=(p.vertices[0:1]), name='block')
    p.Set(vertices=(p.vertices[1:2]), name='force')
    p.Set(edges=(p.edges[3:4], p.edges[5:6]), name='int')
    p.Set(edges=(p.edges[1:2], p.edges[7:8]), name='ext')
    p.Set(edges=(p.edges[2:3], p.edges[6:7]), name='sym')

    # --------------------------------------------------------------------
    # assembly
    # --------------------------------------------------------------------

    a = m.rootAssembly
    a.Instance(dependent=ON, name='base', part=p)

    csys = a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL,
                                    name='Cylindrical',
                                    origin=(0.0, 0.0, 0.0),
                                    point1=(1.0, 0.0, 0.0),
                                    point2=(0.0, 1.0, 0.0))

    # --------------------------------------------------------------------
    # Material
    # --------------------------------------------------------------------

    m.Material(name='hard')
    m.materials['hard'].Density(
        table=((1.04e-09,),))
    m.materials['hard'].Elastic(
        table=((hard_mat, 0.35),))

    m.Material(name='soft')
    m.materials['soft'].Density(
        table=((1.2e-09,),))
    if not hyperelastic:
        m.materials['soft'].Elastic(
            table=((soft_mat, 0.45),))
    else:
        m.materials['soft'].Hyperelastic(materialType=ISOTROPIC, table=((c_coefficient, 0.0),), testData=OFF,
                                         type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

    # --------------------------------------------------------------------
    # Section
    # --------------------------------------------------------------------

    m.HomogeneousShellSection(idealization=NO_IDEALIZATION,
                              integrationRule=SIMPSON,
                              material='hard',
                              name='hard',
                              nodalThicknessField='',
                              numIntPts=5,
                              poissonDefinition=DEFAULT,
                              preIntegrate=OFF,
                              temperature=GRADIENT,
                              thickness=0.1,
                              thicknessField='',
                              thicknessModulus=None,
                              thicknessType=UNIFORM,
                              useDensity=OFF)

    m.HomogeneousShellSection(idealization=NO_IDEALIZATION,
                              integrationRule=SIMPSON,
                              material='soft',
                              name='soft',
                              nodalThicknessField='',
                              numIntPts=5,
                              poissonDefinition=DEFAULT,
                              preIntegrate=OFF,
                              temperature=GRADIENT,
                              thickness=0.1 * rel_depth,
                              thicknessField='',
                              thicknessModulus=None,
                              thicknessType=UNIFORM,
                              useDensity=OFF)

    p.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                        region=p.sets['hard'], sectionName='hard', thicknessAssignment=FROM_SECTION)

    p.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                        region=p.sets['soft'], sectionName='soft', thicknessAssignment=FROM_SECTION)

    # --------------------------------------------------------------------
    # STEP
    # --------------------------------------------------------------------

    # m.ExplicitDynamicsStep(improvedDtMethod=ON, maxIncrement=0.1, name='Analysis', previous='Initial')
    m.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, application=QUASI_STATIC, initialConditions=OFF,
                           initialInc=0.001, minInc=5e-5, maxInc=0.01, maxNumInc=10000, name='Analysis', nlgeom=ON,
                           nohaf=OFF, previous='Initial')
    m.historyOutputRequests['H-Output-1'].setValues(frequency=1)
    m.fieldOutputRequests['F-Output-1'].setValues(frequency=1)
    # m.StaticStep(initialInc=0.001, maxInc=0.01, maxNumInc=10000, minInc=1e-20, name='Analysis', nlgeom=ON, previous='Initial', timePeriod=2.0)
    # m.StaticRiksStep(initialArcInc=0.01, maxArcInc=0.1, maxLPF=1.0, maxNumInc=10000, name='Analysis', nlgeom=ON, previous='Initial')

    # --------------------------------------------------------------------
    # Mesh
    # --------------------------------------------------------------------

    p.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=size)
    p.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges=p.sets['int'].edges, minSizeFactor=0.1, size=size2)
    p.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges=p.sets['ext'].edges, minSizeFactor=0.1, size=size2)
    p.setMeshControls(elemShape=QUAD, regions=p.faces[:], technique=STRUCTURED)
    if quadratic:
        p.setElementType(elemTypes=(ElemType(elemCode=S8R, elemLibrary=STANDARD), ElemType(
            elemCode=STRI65, elemLibrary=STANDARD)), regions=p.sets['all'])
    else:
        p.setElementType(elemTypes=(ElemType(elemCode=S4, elemLibrary=STANDARD), ElemType(
            elemCode=S3, elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=p.sets['all'])

    p.generateMesh()
    a.regenerate()

    # --------------------------------------------------------------------
    # Boundary conditions and loading
    # --------------------------------------------------------------------

    m.DisplacementBC(amplitude=UNSET,
                     createStepName='Analysis',
                     distributionType=UNIFORM,
                     fieldName='',
                     localCsys=None,
                     name='Displacement',
                     region=a.instances['base'].sets['force'],
                     u1=UNSET, u2=UNSET, u3=-0.8 * r_ext, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    m.DisplacementBC(amplitude=UNSET,
                     createStepName='Analysis',
                     distributionType=UNIFORM,
                     fieldName='',
                     localCsys=None,
                     name='Lock',
                     region=a.instances['base'].sets['block'],
                     u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=SET)

    m.DisplacementBC(amplitude=UNSET,
                     createStepName='Analysis',
                     distributionType=UNIFORM,
                     fieldName='',
                     localCsys=a.datums[csys.id],
                     name='Sym',
                     region=a.instances['base'].sets['sym'],
                     u1=UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=UNSET, ur3=SET)

    a.regenerate()

    # --------------------------------------------------------------------
    # job
    # --------------------------------------------------------------------

    j = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model=model_name, modelPrint=OFF,
                multiprocessingMode=DEFAULT, name=project, nodalOutputPrecision=SINGLE,
                numCpus=6, numDomains=6, numGPUs=1, numThreadsPerMpiProcess=1, queue=None,
                resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='',
                waitHours=0, waitMinutes=0)

    mdb.saveAs(cae_file)
    j.submit()
    j.waitForCompletion()
    mdb.saveAs(cae_file)

    odb = session.openOdb(name=odb_name)
    data = session.XYDataFromHistory(name='data', odb=odb, steps=('Analysis',),
                                     outputVariableName='Strain energy: ALLSE for Whole Model')
    del mdb.models[model_name]

    # --------------------------------------------------------------------
    # Deformed shape
    # --------------------------------------------------------------------

    m2 = mdb.Model(modelType=STANDARD_EXPLICIT, name=model_name)
    p3 = m2.PartFromOdb(frame=len(data) - 1, instance='BASE', name='Deformed', odb=session.openOdb(odb_name),
                        shape=DEFORMED, step=0)

    # --------------------------------------------------------------------
    # assembly
    # --------------------------------------------------------------------

    a = m2.rootAssembly
    a.Instance(dependent=ON, name='Deformed', part=p3)
    csys = a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL,
                                    name='Cylindrical',
                                    origin=(0.0, 0.0, 0.0),
                                    point1=(1.0, 0.0, 0.0),
                                    point2=(0.0, 1.0, 0.0))
    # --------------------------------------------------------------------
    # Material
    # --------------------------------------------------------------------

    m2.Material(name='hard')
    m2.materials['hard'].Density(
        table=((1.04e-09,),))
    m2.materials['hard'].Elastic(
        table=((hard_mat, 0.35),))

    m2.Material(name='soft')
    m2.materials['soft'].Density(
        table=((1.2e-09,),))
    if not hyperelastic:
        m2.materials['soft'].Elastic(
            table=((soft_mat, 0.45),))
    else:
        m2.materials['soft'].Hyperelastic(materialType=ISOTROPIC, table=((c_coefficient, 0.0),), testData=OFF,
                                          type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

    # --------------------------------------------------------------------
    # Section
    # --------------------------------------------------------------------

    m2.HomogeneousShellSection(idealization=NO_IDEALIZATION,
                               integrationRule=SIMPSON,
                               material='hard',
                               name='hard',
                               nodalThicknessField='',
                               numIntPts=5,
                               poissonDefinition=DEFAULT,
                               preIntegrate=OFF,
                               temperature=GRADIENT,
                               thickness=0.1,
                               thicknessField='',
                               thicknessModulus=None,
                               thicknessType=UNIFORM,
                               useDensity=OFF)

    m2.HomogeneousShellSection(idealization=NO_IDEALIZATION,
                               integrationRule=SIMPSON,
                               material='soft',
                               name='soft',
                               nodalThicknessField='',
                               numIntPts=5,
                               poissonDefinition=DEFAULT,
                               preIntegrate=OFF,
                               temperature=GRADIENT,
                               thickness=0.1 * rel_depth,
                               thicknessField='',
                               thicknessModulus=None,
                               thicknessType=UNIFORM,
                               useDensity=OFF)

    p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                         region=p3.sets['HARD'], sectionName='hard', thicknessAssignment=FROM_SECTION)

    p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                         region=p3.sets['SOFT'], sectionName='soft', thicknessAssignment=FROM_SECTION)

    # --------------------------------------------------------------------
    # STEP
    # --------------------------------------------------------------------

    # m.ExplicitDynamicsStep(improvedDtMethod=ON, maxIncrement=0.1, name='Analysis', previous='Initial')
    m2.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, application=QUASI_STATIC, initialConditions=OFF,
                            initialInc=0.001, minInc=5e-5, maxInc=0.01, maxNumInc=10000, name='Analysis', nlgeom=ON,
                            nohaf=OFF, previous='Initial')
    m2.historyOutputRequests['H-Output-1'].setValues(frequency=1)
    m2.fieldOutputRequests['F-Output-1'].setValues(frequency=1)
    # m.StaticStep(initialInc=0.001, maxInc=0.01, maxNumInc=10000, minInc=1e-20, name='Analysis', nlgeom=ON, previous='Initial', timePeriod=2.0)
    # m.StaticRiksStep(initialArcInc=0.01, maxArcInc=0.1, maxLPF=1.0, maxNumInc=10000, name='Analysis', nlgeom=ON, previous='Initial')

    # --------------------------------------------------------------------
    # Boundary conditions and loading
    # --------------------------------------------------------------------

    m2.DisplacementBC(amplitude=UNSET,
                      createStepName='Analysis',
                      distributionType=UNIFORM,
                      fieldName='',
                      localCsys=None,
                      name='Displacement',
                      region=a.instances['Deformed'].sets['FORCE'],
                      u1=UNSET, u2=UNSET, u3=1.6 * r_ext, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    m2.DisplacementBC(amplitude=UNSET,
                      createStepName='Analysis',
                      distributionType=UNIFORM,
                      fieldName='',
                      localCsys=None,
                      name='Lock',
                      region=a.instances['Deformed'].sets['BLOCK'],
                      u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=SET)

    m2.DisplacementBC(amplitude=UNSET,
                      createStepName='Analysis',
                      distributionType=UNIFORM,
                      fieldName='',
                      localCsys=a.datums[csys.id],
                      name='Sym',
                      region=a.instances['Deformed'].sets['SYM'],
                      u1=UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=UNSET, ur3=SET)

    a.regenerate()

    # --------------------------------------------------------------------
    # job
    # --------------------------------------------------------------------
    del mdb.jobs['waterbomb']

    j = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model=model_name, modelPrint=OFF,
                multiprocessingMode=DEFAULT, name=project, nodalOutputPrecision=SINGLE,
                numCpus=6, numDomains=6, numGPUs=1, numThreadsPerMpiProcess=1, queue=None,
                resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='',
                waitHours=0, waitMinutes=0)

    mdb.saveAs(cae_file)
    j.submit()
    j.waitForCompletion()
    mdb.saveAs(cae_file)
