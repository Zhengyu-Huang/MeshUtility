import numpy as np
import sys
import FluidGeo

def pair_sort(a,b):
    return (a,b) if a < b else (b,a)



def node_type(n):
    '''
    return node type
    :param n: node id for bumpy wing mesh
    :return: 0:lead edge
             1:top
             2:trialing edge
             3:bottom
    '''
    if(n>=484 and n <550):
        return 0
    elif(n>=140 and n < 484):
        return 3
    elif(n>=115 and n < 140):
        return 2
    else:
        return 1


def read_nodes():
    '''
    Read garbagy file, and get the curve of the birdwing
    :return:
    '''
    mshfile = 'surf.top'
    try:
        fid = open(mshfile, "r")
    except IOError:
        print("File '%s' not found." % mshfile)
        sys.exit()

    nodes = []
    elems = []

    print('Reading mesh ...')

    line = fid.readline()  # should be "Nodes FluidNodes"

    while line:
        line = fid.readline()
        data = line.split()
        if data[0] == 'Elements':
            break;
        nodes.append(list(map(float, data[1:4])))

    while line:
        line = fid.readline()
        if not line:
            break
        data = line.split()
        elems.append(list(map(int, data[2:5])))

    nodes = np.array(nodes, dtype=float)
    elems = np.array(elems, dtype=int) - 1  # fixed top file 1-based index
    fid.close()

    #compute elems nodes
    nNodes = len(nodes) # 0,1...nNodes-1
    elemNodes = set()
    edges = set()
    for ele in elems:
        for i in range(-1,2):
            if ele[i] < nNodes:
                elemNodes.add(ele[i])
                if ele[i+1] < nNodes:
                    edges.add(pair_sort(ele[i],ele[i+1]))


    #output the gmsh file for the 2D bird wing
    elemNodes = list(elemNodes)

    birdShapeNodeNum = len(elemNodes)
    birdShape = np.empty((birdShapeNodeNum, 3))
    for i in range(len(elemNodes)):
        birdShape[i,:] =  nodes[elemNodes[i], :]

    return birdShape

##################################################
#Generate embedded surface
#################################################
def surf_geo(nodes, surf_geo ='embeddedSurf.geo'):
    '''
    :param nodes: a node list, n[0],n[1] ....n[m-1] n[0] is a close curve
    :param surf_geo:
    :return:
    '''
    #0:lead edge,  1:top,   2:trialing edge,   3:bottom
    cl = [0.0005,0.005,0.0005,0.01];
    cl_layer = 0.02;
    layer_n = 50
    file = open(surf_geo, 'w')
    for i in range(4):
        file.write('cl%d = %.15f;\n' %(i,cl[i]))
    file.write('cl_layer = %.15f;\n' %(cl_layer))
    nNodes = len(nodes)
    for i in range(nNodes):

        x,y,z = nodes[i,:]
        nType = node_type(i+1)
        file.write('Point(%d) = {%.15f,   %.15f,   %.15f,   cl%d};\n' %(i+1, x, y, 0.0,nType))

    nLines= nNodes
    for lineId in range(1,nNodes+1):
        file.write('Line(%d) = {%d,   %d};\n' % (lineId, lineId, 1 if lineId == nNodes else lineId+1))

    #extrude line loop
    file.write('surface[] = Extrude {0, 0, %d*'%layer_n)
    file.write('cl_layer')
    file.write('} {Line{')
    file.write(','.join(str(x) for x in range(1,nLines+1)))
    file.write('}; Layers {%d};};\n' %(layer_n))

    file.write('Physical Surface("StickFixedSurface_1") = {')
    file.write(','.join(str(x) for x in range(nLines + 4, nLines + 4*(nLines + 1),4)))
    file.write('};\n')
    file.close()


def domain_geo(birdWing, domain_geo ='domain.geo'):
    '''
    :param birdWing: a node list, n[0],n[1] ....n[m-1] n[0] is a close curve
    :param domain_geo:
    :return:
    '''
    ##################################################
    #Generate fluid domain
    #################################################
    #Background
    x_bg_ball = 0.0
    y_bg_ball = 0.0
    r_bg_ball = 50.0
    cl_bg = 5

    #refinement ball
    x_ref_ball = 1.0
    y_ref_ball = 0.0
    r_ref_ball = 1.2
    cl_ref_ball = 0.01
    d_min_ball = r_ref_ball
    d_max_ball = r_bg_ball


    #refine birdwing
    #0:lead edge,  1:top,   2:trialing edge,   3:bottom
    cl_ref = 0.0005 #[0.0005,0.0005,0.0005,0.01];
    d_min_birdwing = 5*cl_ref
    d_max_birdwing = r_bg_ball



    cl_layer = 0.02;
    layer_n = 50






    domain_geo = 'domain.geo'
    file = open(domain_geo, 'w')
    file.write('cl_ref = %.15f;\n' % (cl_ref))
    file.write('cl_bg = %.15f;\n' % (cl_bg))
    file.write('cl_ref_ball = %.15f;\n' % (cl_ref_ball))
    file.write('cl_layer = %.15f;\n' % (cl_layer))
    file.write('layer_n = %d;\n' % (layer_n))
    file.write('r_ref_ball = %.15f;\n' % (r_ref_ball))


    point_id = 1
    nBirdWing = len(birdWing)
    # point for field 1, bird wing
    for i in range(nBirdWing):
        #disturb = (np.random.uniform(0., 1., 2) - 0.5) * min(cl_capsule, cl_cable) / 2.0  # randomize the points location
        point_id = FluidGeo.writeElem(file, 'Point', point_id, [birdWing[i,0] , birdWing[i,1] , 0.0, 'cl_ref'])



    #point for field 2, ball
    point_id = FluidGeo.writeElem(file, 'Point', point_id, [x_ref_ball, y_ref_ball, 0.0 , 'cl_ref_ball'])




    FluidGeo.writeMeshSize(file, 'Attractor',1, [nBirdWing+1])#ball attractor
    FluidGeo.writeMeshSize(file, 'Attractor', 2, range(1,nBirdWing+1))



    #file, type, field_id, IField, LcMin, LcMax, DistMin, DistMax, Sigmoid
    FluidGeo.writeMeshSize(file, 'Threshold', 3, 1, 'cl_ref_ball', 'cl_bg',           d_min_ball,         d_max_ball,    0) # ball
    FluidGeo.writeMeshSize(file, 'Threshold', 4, 2, 'cl_ref', 'cl_bg',           d_min_birdwing ,        d_max_birdwing, 0) # birdwing
    FluidGeo.writeMeshSize(file, 'Min', 5, [3,4])

    #Background Mesh, the domain is a circle
    FluidGeo.backgroundMesh(file, 'circle', x_bg_ball, y_bg_ball, r_bg_ball, 0.0,  layer_n, point_id)
    file.close()





if __name__ =='__main__':
    nodes = read_nodes()
    surf_geo(nodes)
    domain_geo(nodes)


