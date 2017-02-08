'''
@author: Alister Maguire

A node and edge class to be used within the NewickTree 
class. This node type is specialized for visualization
purposes. 

'''

class Node():
    
    def __init__(self, name="", coords=None):
        self.name   = name
        self.coords = [0.0]*3
        self.left   = None
        self.right  = None
        self.parent = None

        #the following 3 variables are for
        #use within coordinate calculations
        self.degree = 0.0
        self.offset = (0.0, 0.0)
        self.coeff  = (0.0, 0.0)

        self.depth  = 1
        self.e_weights  = [0.0]*3    #let the indices 0, 1, 2 
        self.neighbors  = [None]*3   #map to parent, left, right
                                     #for both weights and neighbors.
        if coords!= None:
            if (type(coords) != list):
                print("ERROR: invalid coords type")
                return 0
            if len(coords) < 3:
                print("ERROR: coods must have 3 elements")
                return 0
            
            for i in range(3):
                self.coords[i] = coords[i]
         
    def get_weights(self):
        return self.e_weights

    def get_coords(self):
        return self.coords
    
    def get_coord(self, i):
        if i > 2 or i < 0:
            print("ERROR: attempting to access coordinate beyond x, y, z")
            return None
        return self.coords[i]

    def get_name(self):
        return self.name

    def get_degree(self):
        return self.degree

    def get_offset(self):
        return self.offset

    def get_coefficient(self):
        return self.coeff

    def get_right(self):
        return self.right

    def get_left(self):
        return self.left
    
    def get_parent(self):
        return self.parent

    def get_neighbors(self):
        return self.neighbors

    def get_depth(self):
        return self.depth

    def set_right(self, right):
        if self.right == None:
            self.degree += 1
        if right == None:
            self.degree -= 1
        else:
            right.depth = self.depth + 1
        self.right = right
        self.neighbors[2] = right
    
    def set_left(self, left):
        if self.left == None:
            self.degree += 1
        if left == None:
            self.degree -= 1
        else:
            left.depth = self.depth + 1
        self.left = left
        self.neighbors[1] = left
    
    def set_parent(self, parent):
        if self.parent == None:
            self.degree += 1
        if parent == None:
            self.degree -= 1
        else:
            #replacing parent => depth remains the same
            parent.depth = self.depth - 1
        self.parent = parent
        self.neighbors[0] = parent

    def set_name(self, name):
        self.name = name
    
    def set_weights(self, weights):
        self.e_weights = weights
    
    def insert_weight(self, w, i):
        self.e_weights[i] = w

    def set_coefficient(self, coeff):
        self.coeff = coeff

    def set_offset(self, offset):
        self.offset = offset

    def set_coord(self, i, val):
        self.coords[i] = val

    def is_leaf(self):
        return True if (self.right == None and self.left == None) else False


class Edge():
    '''
    And edge object that exists as a connection between
    two nodes. 
    '''

    def __init__(self, parent, child):
        self.parent   = parent
        self.child    = child
        self.p_coords = [0.0]*3
        self.c_coords = [0.0]*3
        self.init_coords()

    def init_coords(self):
        p_coords = self.parent.get_coords()
        c_coords = self.child.get_coords() 
        self.p_coords[0] = p_coords[0]
        self.p_coords[1] = p_coords[1]
        self.c_coords[0] = c_coords[0]
        self.c_coords[1] = c_coords[1]

    def get_parent_coords(self):
        return self.p_coords
    
    def get_child_coords(self):
        return self.c_coords
