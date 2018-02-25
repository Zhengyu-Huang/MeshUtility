'''
This is a file, reading parachute structure file and writing parachute embedded file
'''

struFile = 'mesh.in'
import sys
import numpy as np

'''
This is a file, reading parachute structure file and writing parachute embedded file
'''

def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def ReadNodes(file):
    line = file.readline()  # should be "Nodes FluidNodes"
    while line:
        data = line.split()
        if data[0] == 'NODES':
            break
        line = file.readline()
    print('ReadNodes, the first line is ', line)
    nodes = []
    while line:
        line = file.readline()
        data = line.split()
        if data[0][0] == '*' or data[0] == 'NODES':
            continue
        if RepresentsInt(data[0]):
            nodes.append(list(map(float, data[1:4])))
        else:
            break


    return file, line, nodes

def ReadElems(file, line):
    print('ReadElems, the first line is ', line)
    elems = []
    type = -1
    while line:
        line = file.readline()
        data = line.split()
        if not line or data[0][0] == '*':
            continue
        if RepresentsInt(data[0]):
            type = len(data) - 2
            elems.append(list(map(int, data[2:])))

        else:
            break

    return file, line, elems, type


def SplitLines(line_elems):
    '''
    :param line_elems:
    :return:
    '''
    j_old = -1
    line, lines = [], []
    for i,j in line_elems:
        if i != j_old:
            if line:
                lines.append(line)
            line = [i,j]

        else:
            line.append(j)
        j_old = j
    lines.append(line)
    return lines




def LineDressing(line_coord, r, shape):
    '''
    :param line_coord:
    :param r:
    :param shape:
    :param skip:
    :return:
    '''
    n = np.shape(line_coord)[0]
    phantom_coord = np.empty(shape=[2 + n   * shape, 3], dtype=float)
    phantom_tris =  np.empty(shape=[2*shape*n, 3], dtype=float)
    A = line_coord[0, :]
    B = line_coord[-1 , :]
    nx, ny, nz = dir = (B - A) / np.linalg.norm(B - A)  # direction of AB
    e1 = np.array([1.0, 0.0, 0.0], dtype=float)
    e2 = np.array([0.0, 1.0, 0.0], dtype=float)
    e3 = np.array([0.0, 0.0, 1.0], dtype=float)

    theta = 2 * np.pi / shape

    # Rotation matrix on https://en.wikipedia.org/wiki/Rotation_matrix
    R = np.array([[np.cos(theta) + nx * nx * (1 - np.cos(theta)), nx * ny * (1 - np.cos(theta)) - nz * np.sin(theta),
                   nx * nz * (1 - np.cos(theta)) + ny * np.sin(theta)],
                  [ny * nx * (1 - np.cos(theta)) + nz * np.sin(theta), np.cos(theta) + ny * ny * (1 - np.cos(theta)),
                   ny * nz * (1 - np.cos(theta)) - nx * np.sin(theta)],
                  [nz * nx * (1 - np.cos(theta)) - ny * np.sin(theta),
                   nz * ny * (1 - np.cos(theta)) + nx * np.sin(theta), np.cos(theta) + nz * nz * (1 - np.cos(theta))]],
                 dtype=float)

    phantom_dr = np.empty(shape=[shape, 3], dtype=float)

    if(np.linalg.norm(np.cross(dir, e1)) > 0.5):
        phantom_dr[0, :] = r * np.cross(dir, e1)/np.linalg.norm(np.cross(dir, e1))
    elif(np.linalg.norm(np.cross(dir, e2)) > 0.5):
        phantom_dr[0, :] = r * np.cross(dir, e2)/np.linalg.norm(np.cross(dir, e2))
    elif (np.linalg.norm(np.cross(dir, e3)) > 0.5):
        phantom_dr[0, :] = r * np.cross(dir, e3)/np.linalg.norm(np.cross(dir, e3))
    else:
        print('error in LineDressing')


    phantom_coord[0,:],phantom_coord[-1,:] = A,B
    for j in range(1, shape):
        phantom_dr[j, :] = np.dot(R, phantom_dr[j - 1, :])

    for i in range(n):
        for j in range(shape):
            phantom_coord[i * shape + j + 1, :] = line_coord[i, :] + phantom_dr[j, :]


    #Bottom
    for j in range(shape):
        phantom_tris[j , :] = 0,  (j+1)%shape + 1 , j+1

    for i in range(n - 1):
        for j in range(shape):
            phantom_tris[shape + 2 * i * shape + 2 * j, :]     = shape*i + j + 1 , shape*i + (j+1)%shape + 1 , shape*i + (j+1) + shape
            phantom_tris[shape + 2 * i * shape + 2 * j + 1, :] = shape*i + (j+1)%shape + 1 , shape*i + shape + (j+1)%shape + 1 , shape*i + (j+1) + shape



    #Top
    for j in range(shape):
        phantom_tris[-1 - j, :] = shape*n + 1,  shape*n - shape + j + 1, shape*n - shape + (j + 1)%shape + 1

    return phantom_coord, phantom_tris




skip = 1
shape = 4;
r = 0.01;

try:
    struFile = open('./mesh_emb_raw.top', "r")
except IOError:
    print("File '%s' not found." % 'mesh_emb_raw.top')
    sys.exit()

nodes = []
embSurfs = []
bunchLines = []




print('Reading mesh ...')
file, line, nodes = ReadNodes(struFile)
while line:
    '''
    It should read Band Gores, Disk Gores, Gap Lines, Suspension Lines, Vent Lines
    '''
    data = line.split()
    struFile, line, elems, type = ReadElems(struFile, line)
    print('Elem type is ', type)
    if(type == 3):
        embSurfs.append(elems)
    elif(type == 2):
        bunchLines.append(elems)
struFile.close()




print('Processing data ...')
embNodeSet = set()
for elems in embSurfs:
    for i,j,k in elems:
        embNodeSet.add(i)
        embNodeSet.add(j)
        embNodeSet.add(k)
print('Process data ... #Node is ', len(embNodeSet) )
embNodes = list(embNodeSet)
sorted(embNodes)

allLines = []
phantomCoords = [[],[],[]]
phantomTris = [[],[],[]]
for i in range(len(bunchLines)):
    print('line banch ', i)
    lines = SplitLines(bunchLines[i])
    allLines.append(lines)

for i in range(3):
    for line in allLines[i]:
        line_coord = np.empty([len(line) - 2*skip,3], dtype=float)
        for j in range(len(line) - 2*skip):
            line_coord[j,:] = nodes[line[j + skip] - 1]
        phantomCoord, phantomTri = LineDressing(line_coord, r, shape)

        phantomCoords[i].append(phantomCoord)
        phantomTris[i].append(phantomTri)



print('Writing mesh ...')
embFile = open('embeddedSurface.top', 'w')
embFile.write('Nodes nodeset\n')
fabricNodes = len(embNodes)
for i in range(len(embNodes)):
    embFile.write('%d  %.12f  %.12f  %.12f\n' % (i + 1, nodes[i][0], nodes[i][1], nodes[i][2]))

nId = fabricNodes
#write nodes
for i in range(len(phantomCoords)):
    for j in range(len(phantomCoords[i])):
        phantomCoord = phantomCoords[i][j]
        for k in range(len(phantomCoord)):
            nId += 1
            embFile.write('%d  %.12f  %.12f  %.12f\n' % (nId, phantomCoord[k,0], phantomCoord[k,1], phantomCoord[k,2]))



embNodesMap = {}
for i in range(len(embNodes)):
    embNodesMap[embNodes[i]] = i+1
nS = 0;
for nType in range(len(embSurfs)):
    print('nType is ', nType)
    embFile.write('Elements StickMovingSurface_%d using nodeset\n' %(nType+1))
    for i in range(len(embSurfs[nType])):
        nS += 1
        embFile.write('%d  4  %d  %d  %d\n' % (nS,  embNodesMap[embSurfs[nType][i][0]], embNodesMap[embSurfs[nType][i][1]], embNodesMap[embSurfs[nType][i][2]]))

firstNode = fabricNodes + 1
for i in range(len(phantomCoords)):
    nType += 1
    print('nType is ', nType)
    embFile.write('Elements StickMovingSurface_%d using nodeset\n' % (nType+1))
    for j in range(len(phantomTris[i])):

        phantomCoord = phantomCoords[i][j]
        phantomTri = phantomTris[i][j]
        for k in range(len(phantomTri)):
            nS += 1
            embFile.write('%d  4 %d  %d  %d\n' % (nS, firstNode + phantomTri[k,0], firstNode + phantomTri[k,1], firstNode + phantomTri[k,2]))
        firstNode += len(phantomCoord)

embFile.close()



