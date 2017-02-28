#!/usr/bin/python
'''
@author: Alister Maguire

NewickTree() builds a binary tree from a given newick string. 
The tree is specialized to be used for drawing circular, 
phylogenetic trees. Each node within the tree has x, y, z 
coordinates, and these coordinates are computed using a tree
drawing algorithm. 

recognitions:
    Thanks to Christian Bachmaier, Ulrik Brandes, and Barbara Schieper
    for their circular tree drawing algorithm found in their paper titled
    "Drawing Phylogenetic Trees*".

    Thanks to John Conery from the University of Oregon for guidance in 
    handling phylogenetic trees. 

'''

from node import Node, Edge
from counts_map import CountsMap
import argparse
import math

class NewickTree():

    def __init__(self, newick, counts_map):
        self.total_leaves = 0
        self.root         = Node()
        self.counts_map   = counts_map

        #A list of odd characters I've found in newick strings. 
	self.odd_chars    = ['.', ' ', '_', '/', '-', '+', '*']
        self.nodes        = [self.root]
        self.edges        = []
        self.height       = 1
        self.parse_newick(newick)
        self.init_sphere_coordinates()
        self.finalize_coordinates()

    def get_num_samples(self):
        return self.counts_map.get_sample_count()

    def get_root(self):
        return self.root

    def update_depth(self):
        '''
        Update the depth of all nodes. This method
        relies on cascade_depth, which does the actual
        computation. 
        '''
        self.root.depth = 1
        return self.cascade_depth(self.root)

    def cascade_depth(self, parent):
        '''
        Cascade through the nodes and update the depths
        accordingly.
        '''
        if parent.get_left() != None:
            d1 = parent.depth + 1
            parent.get_left().depth = d1
            d2 = self.cascade_depth(parent.get_left())
            return d2 if d2 >= d1 else d1
        if parent.get_right() != None:
            d1 = parent.depth + 1
            parent.get_right().depth = d1
            d2 = self.cascade_depth(parent.get_right())
            return d2 if d2 >= d1 else d1
        

    def preorder(self):
        self.display_preorder(self.root)
    
    def postorder(self):
        self.display_postorder(self.root)

    def inorder(self):
        self.display_inorder(self.root)

    def display_preorder(self, node):
        print(node.get_name())
        if node.get_left() != None:
            self.display_preorder(node.get_left())
        if node.get_right() != None:
            self.display_preorder(node.get_right())

    def display_inorder(self, node):
        if node.get_left() != None:
            self.display_inorder(node.get_left())
        print(node.get_name())
        if node.get_right() != None:
            self.display_inorder(node.get_right())

    def display_postorder(self, node):
        if node.get_left() != None:
            self.display_postorder(node.get_left())
        if node.get_right() != None:
            self.display_postorder(node.get_right())
        print(node.get_name())

    def get_nodes(self):
        return self.nodes
    
    def get_edges(self):
        return self.edges

    def get_sub_str(self, s, i):
        '''
        Get an organism name from the newick string. 
        '''
        sub_str = ''
        counter = 0
        while s[i].isalpha() or s[i] in self.odd_chars:
            sub_str = sub_str + s[i]
            i += 1
            counter += 1

            #I keep finding all sorts of strange characters, 
            #so I've resorted to checking for them now. There's
            #probably a better way of handling this.  
            if counter > 500:
                print("possibly found another invalid character")
                print("char: " + str(s[i]))
                print("exiting...")
                return 0

        #decrement i by 1 because the for loop will
        #re-increment i immediately after. 
        return (sub_str, i-1)

    def parse_newick(self, newick):
        '''
        Parse the newick string, and build our binary
        tree from its structure.  
        '''
        cur_node = self.root #The root is set to an empty node. 
                             #re-rooting takes place at the end. 
        size = len(newick)
        i    = 0

        #Parse the newick string
        while i < size:
            if newick[i] == '(':
                new_left = Node()
                cur_node.set_left(new_left)
                new_left.set_parent(cur_node)
                cur_node = new_left
                self.nodes.append(cur_node)

            elif newick[i] == ')':
                if cur_node.get_parent() != None:
                    cur_node = cur_node.get_parent()
                else:
                    print("ERROR: node missing parent! exiting...")
                    return 0

            elif newick[i] == ',':
                if cur_node.get_parent() != None:
                    cur_node = cur_node.get_parent()
                    new_right = Node()     
                    cur_node.set_right(new_right)
                    new_right.set_parent(cur_node)
                    cur_node = new_right
                    self.nodes.append(cur_node)
                else:
                    print("ERROR: node missing parent! exiting...")
                    return 0

            elif newick[i] == ';':

                #re-root the tree (the initial root is an empty node)
                if self.root.get_left() != None and self.root.get_right() != None:
                    print("ERROR: re-rooting failed; too many children...")
                    return 0
                elif self.root.get_left() != None:
                    self.root = self.root.get_left()
                    self.root.set_parent(None)
                    self.nodes = self.nodes[1:]
                elif self.root.get_right() != None:
                    self.root = self.root.get_right()
                    self.root.set_parent(None)
                    self.nodes = self.nodes[1:]
                else:
                    print("ERROR: re-rooting failed; empty tree..")
                    return 0
                self.height = self.update_depth()
                return 1

            elif newick[i] == ' ':
                print("ERROR: invalid newick character!")
                return 0

            else:
                name, i = self.get_sub_str(newick, i)
                name = name.strip()
                cur_node.set_name(name)
           
            i += 1

    def postorder_circle(self, node, leaves_found):
        '''
        This is the second step in the circular tree algorithm, 
        and it's where most of the computations take place. 
        '''
        
        if node.get_left() != None:
           leaves_found = self.postorder_circle(node.get_left(), leaves_found)
        if node.get_right() != None:
           leaves_found = self.postorder_circle(node.get_right(), leaves_found)

        if node.is_leaf() or (self.root.get_name() 
           == node.get_name() and node.get_degree() == 1):
            #put leaf nodes on the unit circle
            node.set_coefficient((0.0, 0.0))
            const  = (2.0*math.pi*leaves_found)/self.total_leaves 
            offset = (math.cos(const), math.sin(const))
            node.set_offset(offset)
            leaves_found += 1

        else:
            #map coordinates for internal nodes
            e_total   = 0.0
            e_weights = [0.0]*3
            neighbors = node.get_neighbors()
           
            for i in range(3):
                if neighbors[i] != None:
                    if (node.get_name() == self.root.get_name()
                       or i == 0):
                        e_weights[i] = 1.0/node.get_depth()
                    else:
                        denom = float(node.get_depth())*float((node.get_degree()-1))
                        e_weights[i] = 1.0/denom
                    e_total += e_weights[i]

            node.set_weights(e_weights)

            x_t  = 0.0
            x_t2 = 0.0
            y_t  = 0.0
            y_t2 = 0.0

            for i in range(1, 3):            
                if neighbors[i] != None:
                    x_t  += (e_weights[i]/e_total)*(neighbors[i].get_coefficient()[0])
                    x_t2 += (e_weights[i]/e_total)*(neighbors[i].get_offset()[0])
                    y_t  += (e_weights[i]/e_total)*(neighbors[i].get_coefficient()[1])
                    y_t2 += (e_weights[i]/e_total)*(neighbors[i].get_offset()[1])
                                                                              
            if node.get_name() != self.root.get_name():
                x_co = e_weights[0]/(e_total*(1.0-x_t))
                y_co = e_weights[0]/(e_total*(1.0-y_t))
                node.set_coefficient( (x_co, y_co) )

            node.set_offset((x_t2/(1.0-x_t), (y_t2/(1.0-y_t))))
        return leaves_found
        
    
    def preorder_circle(self, node):
        '''
        The final step in the circular 
        tree algorithm. 
        '''
        if node.get_name() == self.root.get_name():
            #Set the roots coordinates
            node.set_coord(0, node.get_offset()[0]) 
            node.set_coord(1, node.get_offset()[1]) 
        else:
            #Set coordinates for all other nodes and leaves
            x = ((node.get_coefficient()[0]*node.get_parent().get_coord(0))
                 + node.get_offset()[0])
            y = ((node.get_coefficient()[1]*node.get_parent().get_coord(1))
                 + node.get_offset()[1])
            node.set_coord(0, x) 
            node.set_coord(1, y) 

        if node.get_left()!= None:
            self.preorder_circle(node.get_left())
        if node.get_right()!= None:
            self.preorder_circle(node.get_right())

    def init_sphere_coordinates(self):
        '''
        This is the circular tree algorithm noted in
        the header. x, y coordinates are calculated for
        each node in the tree. Leaves are place on the unit
        circle, and inner nodes are placed within the unit
        circle. The root is generally skewed to one side, though
        there is a re-rooting option (different from the simple 
        re-rooting I am performing). 
        '''
        leaves_found = 0
        for n in self.nodes:
            if n.get_degree() == 1:
                self.total_leaves += 1 

        self.postorder_circle(self.root, leaves_found)
        self.preorder_circle(self.root)

    def finalize_coordinates(self):
        '''
        The initial coordinates end up on the unit circle. 
        If the number of leaves > 5, I scale the coordinates
        by 10% of the total number of leaves. 
        Also, this is where the edge list is computed, as the
        edge coordinates must be computed after the node 
        coordinates are finalized. 
        '''
        size = len(self.nodes)
        for i in range(size):
            cur_node = self.nodes[i]
            if cur_node.is_leaf():
                counts = self.counts_map.get_counts(cur_node.get_name())
                if counts == None:
                    print('ERROR: leaf name not in counts_map!')
                    print('Aborting tree creation')
                    return 0
                self.nodes[i].set_counts(counts)

            x = cur_node.get_coord(0)
            y = cur_node.get_coord(1)
            if self.total_leaves > 5:
                scale = self.total_leaves*.1 
                self.nodes[i].set_coord(0, x*scale)
                self.nodes[i].set_coord(1, y*scale)
            else:
                self.nodes[i].set_coord(0, x)
                self.nodes[i].set_coord(1, y)
            
            #create edge list
            if self.nodes[i].get_parent() != None:
                e = Edge(self.nodes[i].get_parent(), self.nodes[i])
                self.edges.append(e)

         
                
            
if __name__ == '__main__':
    '''
    Command line testing. 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('newick_file', type=str)
    parser.add_argument('condensed_counts_file', type=str)
    args  = parser.parse_args()
    n_f   = open(args.newick_file, "r")

    n_str = n_f.readlines()[0]
    c_map = CountsMap(condensed_counts_file)
    tree  = NewickTree(n_str, c_map)
    tree.preorder()

