'''
This is a file, reading parachute structure file and writing parachute embedded file
'''

struFile = 'mesh.in'
import sys

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






try:
    struFile = open('./mesh_emb_raw.top', "r")
except IOError:
    print("File '%s' not found." % 'mesh_Structural.in.15')
    sys.exit()

nodes = []
embSurfs = []
embSurfTypes= []




print('Reading mesh ...')
file, line, nodes = ReadNodes(struFile)
while line:
    data = line.split()
    struFile, line, elems, type = ReadElems(struFile, line)
    print('Elem type is ', type)
    if(type == 3):
        embSurfs.append(elems)
        embSurfTypes.append(type)
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



print('Writing mesh ...')
embFile = open('embeddedSurface.top', 'w')
embFile.write('Nodes nodeset\n')
for i in range(len(embNodes)):
    embFile.write('%d  %.12f  %.12f  %.12f\n' % (i + 1, nodes[i][0], nodes[i][1], nodes[i][2]))

embNodesMap = {}
for i in range(len(embNodes)):
    embNodesMap[embNodes[i]] = i+1

nS = 0;
for nType in range(len(embSurfTypes)):
    embFile.write('Elements StickMovingSurface_%d using nodeset\n' %(nType+1))
    for i in range(len(embSurfs[nType])):
        nS += 1
        embFile.write('%d  4  %d  %d  %d\n' % (nS,  embNodesMap[embSurfs[nType][i][0]], embNodesMap[embSurfs[nType][i][1]], embNodesMap[embSurfs[nType][i][2]]))

embFile.close()



