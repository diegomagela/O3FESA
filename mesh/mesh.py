# %%

import gmsh
import sys
import numpy

# Initialize Gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("plate")

# Geometry definition

L = 11.5  # plate length

# Points (counterclockwise)
point_1 = gmsh.model.geo.addPoint(0, 0, 0)
point_2 = gmsh.model.geo.addPoint(L, 0, 0)
point_3 = gmsh.model.geo.addPoint(L, L, 0)
point_4 = gmsh.model.geo.addPoint(0, L, 0)

# Lines (counterclockwise)
bottom = gmsh.model.geo.addLine(1, 2)
left = gmsh.model.geo.addLine(2, 3)
top = gmsh.model.geo.addLine(3, 4)
right = gmsh.model.geo.addLine(4, 1)

# Curve loop and plane surface
curveLoop = gmsh.model.geo.addCurveLoop([bottom, left, top, right])
planeSurface = gmsh.model.geo.addPlaneSurface([curveLoop])

# Mesh

# Set transfinite lines
# nNodesPerSide = int(sys.argv[1])
nNodesPerSide = 3
transLineBottom = gmsh.model.geo.mesh.setTransfiniteCurve(bottom, nNodesPerSide)
transLineLeft = gmsh.model.geo.mesh.setTransfiniteCurve(left, nNodesPerSide)
transLineTop = gmsh.model.geo.mesh.setTransfiniteCurve(top, nNodesPerSide)
transLineRight = gmsh.model.geo.mesh.setTransfiniteCurve(right, nNodesPerSide)

# Set transfinite surface
transSurface = gmsh.model.geo.mesh.setTransfiniteSurface(
    1, "Left", [point_1, point_2, point_3, point_4]
)

# To create rectangles instead of triangles
gmsh.model.geo.mesh.setRecombine(2, 1)

# Synchronize data
gmsh.model.geo.synchronize()

# Physical groups
boundaries = gmsh.model.addPhysicalGroup(
    1, [top, left, bottom, right], name="Boundaries"
)
leftRight = gmsh.model.addPhysicalGroup(1, [left, right], name="leftRight")
bottomTop = gmsh.model.addPhysicalGroup(1, [bottom, top], name="bottomTop")
surface = gmsh.model.addPhysicalGroup(2, [planeSurface], name="Surface")

# Generate mesh

meshOrder = 1
gmsh.model.mesh.setOrder(meshOrder)
gmsh.model.mesh.generate(2)
meshOrder = 2
gmsh.model.mesh.setOrder(meshOrder)

# Export input file
file = open("mesh.inp.in", "w")  # open file and write input file

# *NODE
# TAG X Y Z

file.write("*NODE" + "\n")

nodeTags = gmsh.model.mesh.getNodes()[0]
nNodes = len(nodeTags)
nodesCoords = numpy.array(gmsh.model.mesh.getNodes(-1, -1)[1])
nodesCoords = nodesCoords.reshape((nNodes, 3))

for i in range(nNodes):
    line = []
    line = [nodeTags[i]]
    line.extend(nodesCoords[i, :].tolist())
    file.write(", ".join(map(str, line)))
    file.write("\n")

# *ELEMENT
# TAG NODE_1 NODE_2 ... NODE_N

file.write("*ELEMENT, TYPE=Q9, ELSET=EALL" + '\n')

firstElementTag = gmsh.model.mesh.getElements(2, -1)[1][0][0]
elementTags = gmsh.model.mesh.getElements(2, -1)[1][0] - firstElementTag + 1
nElements = len(elementTags)
elementsNodes = numpy.array(gmsh.model.mesh.getElements(2, -1)[2])
elementsNodes = elementsNodes.reshape((nElements, 9))
elementsNodes = elementsNodes

for i in range(nElements):
    line = []
    line = [elementTags[i]]
    line.extend(elementsNodes[i, :].tolist())
    file.write(", ".join(map(str, line)))
    file.write("\n")

# *BOUNDARY
# NODE_TAG U_X U_Y U_Z ROT_X ROT_Y ROT_Z VALUE
# Named constraints
# ENCASTRE  Constraint on all displacements and rotations at a node.
# PINNED    Constraint on all translational degrees of freedom.
# XSYMM     Symmetry constraint about a plane of constant x coordinate.
# YSYMM     Symmetry constraint about a plane of constant y coordinate.
# ZSYMM     Symmetry constraint about a plane of constant z coordinate.
# XASYMM    Antisymmetry constraint about a plane of constant x coordinate.
# YASYMM    Antisymmetry constraint about a plane of constant y coordinate.
# ZASYMM    Antisymmetry constraint about a plane of constant z coordinate.

file.write("*BOUNDARY" + '\n')

clampedNodes = gmsh.model.mesh.getNodesForPhysicalGroup(1, boundaries)[0]
nBoundaryNodes = len(clampedNodes)

for i in clampedNodes:
    line = [i, "ENCASTRE"]
    file.write(", ".join(map(str, line)))
    file.write("\n")

# *CLOAD
# NODE_TAG DOF VALUE
file.write("*CLOAD" + '\n')

cLoadNodes = [17]

wLoad = 1

for i in cLoadNodes:
    line = [i, 3, wLoad]
    file.write(", ".join(map(str, line)))
    file.write("\n")

# *DLOAD
# ELEMENT_TAG PRESSURE_VALUE
file.write("*DLOAD" + '\n')


tagPhysicalGroup = gmsh.model.getEntitiesForPhysicalGroup(2, 4)
elementsPhysicalGroup = gmsh.model.mesh.getElements(2, tagPhysicalGroup[0])
firstElementPhysicalGroupTag = elementsPhysicalGroup[1][0][0]
elementTagsPhysicalGroup = (
    gmsh.model.mesh.getElements(2, -1)[1][0] - firstElementPhysicalGroupTag + 1
)
nElementsPhysicalGroup = len(elementTagsPhysicalGroup)

pressure = 1

for i in elementTagsPhysicalGroup:
    line = [i, pressure]
    file.write(", ".join(map(str, line)))
    file.write("\n")

# *SHELL SECTION, ELSET=ELSET_NAME, COMPOSITE
# LAYER THICKNESS, NUMBER OF INTEGRATION POINTS, MATERIAL, ORIENTATION
layerThickness = 0.1
nIntPts = 0
materialName = "CFRP"
angleVec = [0, 90]
file.write("*SHELL SECTION, ELSET=EALL, COMPOSITE" + '\n')

for i in angleVec:
    line = [layerThickness, nIntPts, materialName, i]
    file.write(", ".join(map(str, line)))
    file.write("\n")

# *MATERIAL, NAME=CFRP
# *ELASTIC, TYPE=LAMINA
# E_1, E_2, NU_12, G_12, G_13, G_23, TEMPERATURE
E_1 = 1E9
E_2 = 0.5 * E_1
NU_12 = 0.3
G_12 = 0.25 * E_1
G_13 = 0.1 * E_1
G_23 = G_13
density = 1000


file.write("*MATERIAL, NAME=CFRP" + '\n')
file.write("*ELASTIC, TYPE=LAMINA" + '\n')
line = [E_1, E_2, NU_12, G_12, G_13, G_23]
file.write(", ".join(map(str, line)))
file.write("\n")
file.write("*DENSITY" + '\n')
file.write(str(density) + '\n')

# *MATERIAL, NAME=ISO
# *ELASTIC, TYPE=ISOTROPIC
# E, NU
E = 1E9
NU = 0.3
density = 1000

file.write("*MATERIAL, NAME=ISO" + '\n')
file.write("*ELASTIC, TYPE=ISOTROPIC" + '\n')
line = [E, NU]
file.write(", ".join(map(str, line)))
file.write("\n")
file.write("*DENSITY" + '\n')
file.write(str(density) + '\n')


file.close()


# gmsh.fltk.run()
# %%
