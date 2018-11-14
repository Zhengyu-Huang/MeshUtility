#################################
# StructureFile
# fixed displacement point
# fixed rotation point
#
# mesh_Structural.top.quad
# NODES
#
#
# mesh_Structural.surfacetop.quad
#################################
import sys
import numpy as np
class Elem:
    def __init__(self, id, nodes, att = 0, eframe = None):
        self.id = id
        self.nodes = nodes
        self.att = att
        self.eframe = eframe

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
        if data[0] == 'NODES' or data[0] == 'Nodes':
            break
        line = file.readline()
    print('ReadNodes, the first line is ', line)
    nodes = []
    while line:
        line = file.readline()
        data = line.split()
        if data[0][0] == '*' or data[0] == 'NODES' or data[0] == 'Nodes':
            continue
        if RepresentsInt(data[0]):
            nodes.append(list(map(float, data[1:4])))
        else:
            break

    print("ReadNodes reads ", len(nodes), " nodes")
    return file, line, nodes
def pair(a,b):
    return (min(a,b), max(a,b))

def ReadElems(file, line):
    print('\n\n*ReadElems, the first line is ', line)
    elems = []
    type = -1
    name = line.split()[1]
    while line:
        line = file.readline()
        data = line.split()
        if not line or data[0][0] == '*':
            continue
        if RepresentsInt(data[0]):
            type = len(data) - 2
            elems.append(Elem(int(data[0]), list(map(int, data[2:]))))

        else:
            break
    if line.split()[0]  == 'ATTRIBUTES':
        ind = 0
        while line:
            line = file.readline()
            data = line.split()
            if not line or data[0][0] == '*':
                continue
            if RepresentsInt(data[0]):
                elems[ind].att = int(data[1])
                assert(elems[ind].id == int(data[0]))
                ind += 1

            else:
                break
    if line.split()[0]  == 'EFRAMES':
        ind = 0
        while line:
            line = file.readline()
            data = line.split()
            if not line or data[0][0] == '*':
                continue
            if RepresentsInt(data[0]):
                elems[ind].eframe = list(map(float, data[1:]))
                assert(elems[ind].id == int(data[0]))
                ind += 1

            else:
                break
    print('ReadElems reads ', len(elems), ' ', name, ' elems')
    return file, line, elems, type, name


class Mesh:
    '''
    #  nodes set
    #  topology set
    '''
    def __init__(self):
        self.nodes = None
        self.ele_set = []
        self.ele_set_info = []

    def read_stru(self, file_name):
        '''
        :param inputStru:
        :param beamPars: parameters to handle phantom surface, skip beamPars[0] beams at all ends of these lines, the shape of
        the cross section beamPars[1], 4 for square, radius of the cross section beamPars[2]
        :return: nodes: a list of 3 double,
                 elems: a list of 3 int
        '''

        try:
            stru_file = open(file_name, "r")
        except IOError:
            print("File '%s' not found." % file_name)
            sys.exit()


        print('Reading Structure mesh ...')
        file, line, self.nodes = ReadNodes(stru_file)
        while line:
            '''
            It should read (0)Band Gores, Disk Gores, (1)Gap Lines, (2)Suspension Lines, (3)Vent Lines
            '''
            stru_file, line, elems, type, name = ReadElems(stru_file, line)
            self.ele_set.append(elems)
            self.ele_set_info.append([name, type])
        stru_file.close()

    def refine(self, refine_all_beams_or_not = True):
        '''
        add nodes, append to nodes
        change quad
        change beams in quad
        change other beams
        :return:
        '''


        nodes = self.nodes
        ele_set = self.ele_set
        ele_set_info = self.ele_set_info

        n_n = len(nodes) # save number of nodes
        n_es = len(ele_set)

        edge_to_center_node = {}
        new_ele_set = []

        ################   Update canopy
        for i_es in range(n_es):
            ele = ele_set[i_es]
            ele_info = ele_set_info[i_es]
            new_ele = []
            if ele_info[1] == 2:
                print('beam')
                #beam
            elif ele_info[1] == 4:
                print('quad')
                n_e = len(ele)
                for i_e in range(n_e):
                    id = ele[i_e].id
                    att = ele[i_e].att
                    eframe = ele[i_e].eframe
                    ele_nodes = ele[i_e].nodes
                    new_nodes = [0,0,0,0,0]
                    # step 1: add a new node at the center
                    new_nodes[0] = n_n + 1
                    n_n += 1

                    # update node coordinate
                    # todo
                    new_node_coord = [(nodes[ele_nodes[0] - 1][0] + nodes[ele_nodes[1] - 1][0] + nodes[ele_nodes[2] - 1][0] + nodes[ele_nodes[3] - 1][0]) / 4.0,
                                      (nodes[ele_nodes[0] - 1][1] + nodes[ele_nodes[1] - 1][1] + nodes[ele_nodes[2] - 1][1] + nodes[ele_nodes[3] - 1][1]) / 4.0,
                                      (nodes[ele_nodes[0] - 1][2] + nodes[ele_nodes[1] - 1][2] + nodes[ele_nodes[2] - 1][2] + nodes[ele_nodes[3] - 1][2]) / 4.0]

                    nodes.append(new_node_coord)

                    for i_n in range(ele_info[1]):
                        if pair(ele_nodes[i_n - 1], ele_nodes[i_n]) in edge_to_center_node:
                            #the node is exists
                            new_nodes[i_n + 1] = edge_to_center_node[pair(ele_nodes[i_n - 1], ele_nodes[i_n])]
                        else:
                            new_nodes[i_n + 1] = n_n + 1
                            # update map
                            edge_to_center_node[pair(ele_nodes[i_n - 1], ele_nodes[i_n])] = n_n + 1
                            # update node coordinate
                            #todo
                            new_node_coord = [(nodes[ele_nodes[i_n - 1] - 1][0] + nodes[ele_nodes[i_n] - 1][0]) / 2.0,
                                              (nodes[ele_nodes[i_n - 1] - 1][1] + nodes[ele_nodes[i_n] - 1][1]) / 2.0,
                                              (nodes[ele_nodes[i_n - 1] - 1][2] + nodes[ele_nodes[i_n] - 1][2]) / 2.0]

                            nodes.append(new_node_coord)

                            n_n += 1

                    new_ele.append(Elem(id, [ele_nodes[0], new_nodes[2], new_nodes[0], new_nodes[1]], att, eframe))
                    new_ele.append(Elem(id, [ele_nodes[1], new_nodes[3], new_nodes[0], new_nodes[2]], att, eframe))
                    new_ele.append(Elem(id, [ele_nodes[2], new_nodes[4], new_nodes[0], new_nodes[3]], att, eframe))
                    new_ele.append(Elem(id, [ele_nodes[3], new_nodes[1], new_nodes[0], new_nodes[4]], att, eframe))

            new_ele_set.append(new_ele)

        ################Update beams on the canopy
        for i_es in range(n_es):
            ele = ele_set[i_es]
            ele_info = ele_set_info[i_es]
            new_ele = new_ele_set[i_es]
            if ele_info[1] == 2:
                # beam
                n_e = len(ele)
                for i_e in range(n_e):
                    id = ele[i_e].id
                    att = ele[i_e].att
                    eframe = ele[i_e].eframe
                    ele_nodes = ele[i_e].nodes

                    if pair(ele_nodes[0], ele_nodes[1]) in edge_to_center_node:
                        # the node is exists
                        new_nodes = edge_to_center_node[pair(ele_nodes[0], ele_nodes[1])]
                        # update element
                        # todo attributes and eframe
                        new_ele.append(Elem(id, [ele_nodes[0], new_nodes], att, eframe))
                        new_ele.append(Elem(id, [new_nodes, ele_nodes[1]], att, eframe))


        ################Update beams not on the canopy

        for i_es in range(n_es):
            ele = ele_set[i_es]
            ele_info = ele_set_info[i_es]
            new_ele = new_ele_set[i_es]
            if ele_info[1] == 2:
                # beam
                n_e = len(ele)
                for i_e in range(n_e):
                    id = ele[i_e].id
                    att = ele[i_e].att
                    eframe = ele[i_e].eframe
                    ele_nodes = ele[i_e].nodes

                    if pair(ele_nodes[0], ele_nodes[1]) in edge_to_center_node:
                        # the node is exists
                        continue
                    else:
                        if refine_all_beams_or_not:
                            new_nodes = n_n + 1
                            # update map
                            edge_to_center_node[pair(ele_nodes[0], ele_nodes[1])] = n_n + 1
                            # update node coordinate
                            # todo
                            new_node_coord = [(nodes[ele_nodes[0] - 1][0] + nodes[ele_nodes[1] - 1][0])/2.0,
                                              (nodes[ele_nodes[0] - 1][1] + nodes[ele_nodes[1] - 1][1])/2.0,
                                              (nodes[ele_nodes[0] - 1][2] + nodes[ele_nodes[1] - 1][2])/2.0]
                            nodes.append(new_node_coord)
                            n_n += 1

                            new_ele.append(Elem(id, [ele_nodes[0], new_nodes], att, eframe))
                            new_ele.append(Elem(id, [new_nodes, ele_nodes[1]], att, eframe))
                        else :# do not refine these beams
                            new_ele.append(Elem(id, [ele_nodes[0], ele_nodes[1]], att, eframe))


        self.ele_set = new_ele_set

    def write_stru(self, stru_file_name, surf_file_name, thickness = 2.e-3):
        print('Writing mesh ...')
        stru_file = open(stru_file_name, 'w')
        surf_file = open(surf_file_name, 'w')
        stru_file.write('NODES\n')

        # Step1.1 write nodes
        nodes = self.nodes
        n_n = len(nodes)
        for i in range(n_n):
            stru_file.write('%d  %.16E  %.16E  %.16E\n' % (
            i + 1, nodes[i][0], nodes[i][1], nodes[i][2]))

        # Step1.2 write TOPOLOGY
        ele_set = self.ele_set
        ele_set_info = self.ele_set_info
        n_es = len(ele_set)
        stru_ele_start_id = 1
        surf_ele_start_id = 1
        for i in range(n_es):
            ele_new = ele_set[i]
            ele_info = ele_set_info[i]
            stru_file.write('*  %s\n' % ele_info[0])
            stru_file.write('TOPOLOGY\n')
            n_e = len(ele_new)
            type = fem_type = ele_info[1]
            if type == 2:
                fem_type = 6
            elif type == 4:
                fem_type = 16

            for j in range(n_e):
                if (type == 2):
                    stru_file.write('%d  %d  %d  %d\n' % (stru_ele_start_id + j, fem_type, ele_new[j].nodes[0], ele_new[j].nodes[1]))
                if (type == 4):
                    stru_file.write(
                        '%d  %d  %d  %d  %d  %d\n' % (stru_ele_start_id + j, fem_type, ele_new[j].nodes[0], ele_new[j].nodes[1], ele_new[j].nodes[2], ele_new[j].nodes[3]))
            stru_file.write('ATTRIBUTES\n')
            for j in range(n_e):
                stru_file.write('%d  %d\n' % (stru_ele_start_id + j, ele_new[j].att))

            if(type == 2): #beam elements need EFRAMES
                stru_file.write('EFRAMES\n')
                for j in range(n_e):
                    stru_file.write('%d  %.16E  %.16E  %.16E  %.16E  %.16E  %.16E  %.16E  %.16E  %.16E\n' % (stru_ele_start_id + j,
                                                                                                             ele_new[j].eframe[0], ele_new[j].eframe[1], ele_new[j].eframe[2],
                                                                                                             ele_new[j].eframe[3], ele_new[j].eframe[4], ele_new[j].eframe[5],
                                                                                                             ele_new[j].eframe[6], ele_new[j].eframe[7], ele_new[j].eframe[8]))

            stru_ele_start_id += n_e
            ############Write Sufrace top
            if (type == 4):
                name = ele_info[0]
                surf_id = 1 if name == 'Disk_Gores' else 2
                surf_file.write('SURFACETOPO %d SURFACE_THICKNESS %.16E\n' %(surf_id, thickness))

                for j in range(n_e):
                    surf_file.write(
                        '%d  %d  %d  %d  %d  %d\n' % (
                         j + surf_ele_start_id, 1, ele_new[j].nodes[0], ele_new[j].nodes[1], ele_new[j].nodes[2], ele_new[j].nodes[3]))
                surf_ele_start_id += n_e

        stru_file.close()
        surf_file.close()

    def visualize_disp(self):
        '''
        update node coordinate to include the displacement
        :return:
        '''
        nodes = self.nodes
        node_disp = self.node_disp
        nn = len(nodes)
        for i_n in range(nn):
            nodes[i_n] = nodes[i_n][0] + node_disp[i_n][0], nodes[i_n][1] + node_disp[i_n][1], nodes[i_n][2] + node_disp[i_n][2]

    def folding(self, n):
        '''
        First fold the Disk, then the band, and the suspension lines
        n : number of
        :return: disp

        '''
        nodes = self.nodes
        ele_set = self.ele_set
        ele_set_info = self.ele_set_info
        node_disp = np.zeros([len(nodes),3])
        #parachute parameters
        r_d, R_d, h_d =  0.788, 7.7235, 39.2198 #Disk inner radius, outer radius and z
        R_b, ht_b, hb_b = 7.804, 38.3158, 35.7358 #Band radius, band top z and band bottom z
        h0 = 0

        ############ piece map
        theta = np.pi / n
        # the first parameter is cosa in [0, 1]
        cosa = 0.1
        sina = -np.sqrt(1 - cosa * cosa)
        cosb = (sina * cosa * np.tan(theta) - sina / np.cos(theta)) / (1 + sina * sina * np.tan(theta) * np.tan(theta))
        sinb = -cosa + cosb * sina * np.tan(theta)
        print('cosb^2 + sinb^2 = ', cosb * cosb + sinb * sinb)

        disk_rot0 = np.array([[cosa, -cosb * sina, -sinb * sina],
                         [0, cosa - cosb * sina * np.tan(theta), cosb],
                         [sina, cosb * cosa, sinb * cosa]])
        # the second parameter is b3
        b3 = 50.0
        disk_disp = np.array([0, 0, b3])
        print('orthogonal matrix test ', np.dot(disk_rot0.T, disk_rot0))

        # The rigid motion of the first half gore is rot0*x + disp0
        x1 = np.array([r_d, 0, 0])
        x2 = np.array([R_d, 0, 0])
        x3 = np.array([r_d*np.cos(theta), r_d*np.sin(theta), 0])
        x4 = np.array([R_d*np.cos(theta), R_d*np.sin(theta), 0])

        y1 = np.dot(disk_rot0, x1) + disk_disp
        print('y1 is ', y1)
        y2 = np.dot(disk_rot0, x2) + disk_disp
        print('y1 is ', y2)
        y3 = np.dot(disk_rot0, x3) + disk_disp
        print('y1 is ', y3)
        y4 = np.dot(disk_rot0, x4) + disk_disp
        print('y1 is ', y4)

        disk_rot_matrices = []

        ######################
        for gore_id in range(n):
            #Handle the piece of theta = [2pi/n * gore_id, 2pi/n * (gore_id+1)]
            #The piece will be fold at the center

            ##
            # The the map is
            # rotn = R_z(theta*n) * rot0 * R_z(theta*n)^{-1}
            # bi = b0
            ##
            R_z = np.array([[np.cos(2*gore_id*theta), -np.sin(2*gore_id*theta), 0.],
                            [np.sin(2*gore_id*theta), np.cos(2*gore_id*theta), 0.],
                            [0.,0.,1.]])
            disk_rotn = np.dot(R_z, np.dot(disk_rot0, R_z.T))

            Ref_z = np.array([[np.cos((4*gore_id+2)*theta), np.sin((4*gore_id+2)*theta), 0.],
                              [np.sin((4*gore_id+2)*theta), -np.cos((4*gore_id+2)*theta), 0.],
                              [0., 0., 1.]])
            disk_rotn_2 = np.dot(Ref_z, np.dot(disk_rotn, Ref_z.T))
            disk_rot_matrices.append(disk_rotn)
            disk_rot_matrices.append(disk_rotn_2)

        #For the piece  [2pi/n * gore_id, 2pi/n * (gore_id + 0.5)]
        #map the point (r, theta) to

        n_n = len(nodes)  # save number of nodes
        n_es = len(ele_set)

        ################   Update canopy
        for i_es in range(n_es):
            ele = ele_set[i_es]
            ele_info = ele_set_info[i_es]
            if ele_info[1] == 4 and ele_info[0] == 'Disk_Gores':
                n_e = len(ele)
                for i_e in range(n_e):
                    ele_nodes = ele[i_e].nodes
                    for i_n in ele_nodes:
                        xx = nodes[i_n - 1]
                        angle_x = np.arctan2(xx[1], xx[0]) + 2*np.pi#[-pi, 2pi]
                        gore_id = int(angle_x/theta)%(2*n)


                        disk_rot = disk_rot_matrices[gore_id]
                        xx_shift, disp_shift = np.array([xx[0], xx[1], 0.0]), np.array([0, 0, xx[2]])
                        new_xx = np.dot(disk_rot, xx_shift) + disk_disp + disp_shift

                        node_disp[i_n - 1,:] =  new_xx[0] - xx[0], new_xx[1] - xx[1], new_xx[2] - xx[2]




        self.node_disp = node_disp






if __name__ == '__main__':
    mesh = Mesh()

    mesh.read_stru('Parachute_Quad_Init/mesh_Structural.top.quad')
    mesh.folding(8)
    # mesh.refine()
    mesh.visualize_disp()
    mesh.write_stru('mesh_Structural.top.quad', 'mesh_Structural.surfacetop.quad')







