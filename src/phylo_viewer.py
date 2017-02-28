#!/usr/bin/python
'''
@author: Alister Maguire

The first steps towards creating a 3D phyolgenetic tree
viewer that is to be used for interpreting time sliced 
trees. 
In particular, we are aiming for a means of visualizing 
the results of microbiota studies, where samples of the 
micriobiome populations are measured at different times. 
'''
from node import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGL.GL import shaders
import sys
import math
from newick_tree import NewickTree
from counts_map import CountsMap
import numpy
import argparse
import ctypes

MAX_LAYERS = 30
SPACING    = 3

class TreeViewer():
    '''
       A 3d pyholgenetic tree viewer (under construction) 
    '''
    def __init__(self, nodes, edges, leaf_count, sample_count, layers):
    
        self.layer_count = layers
        self.num_samples = sample_count
        self.leaf_count  = leaf_count
        self.nodes       = nodes
        self.edges       = edges
        self.rot_y_left  = 0
        self.rot_y_right = 0
        self.rot_x_up    = 0
        self.rot_x_down  = 0
        self.rot_z_left  = 0
        self.rot_z_right = 0
        self.y_deg       = 0
        self.x_deg       = 0
        self.z_deg       = 0
        self.zoom_in     = 0
        self.zoom_out    = 0
        self.zoom_val    = 0
        self.start_zoom  = 10

        #Create the circle geometry
        start_z = -1*SPACING*(self.layer_count/2) #(self.num_samples/100)*-5
        raw_points = []
        cur_z      = start_z
        self.cylinders  = []
        self.num_circles = 0
        self.num_edges   = 0

        colors = []

        leaf_color = [0.2, 0.8, 0.3, 1.0]
        #Prototype for leaf interpolation. 
        #Eventually, I'll want to send the cylinder geometry to the
        #gpu for direct rendering. For now, though, I'm just using gluCylinder. 
        for node in self.nodes:
            x, y, z = node.get_coords()
            if node.is_leaf():
                samples = node.get_counts_list() 
                i = 0
                 
                #if we only have one layer, just create spheres. 
                if self.layer_count == 1:
                        self.num_circles += 1
                        raw_points.extend(self.create_circle(20*samples[i], x, y, start_z + i*SPACING))
                        colors.extend(leaf_color*360)

                #if we have multiple layers, look for cylinders.  
                
                while i < (self.layer_count-1):
                    if samples[i] > 0 and samples[i+1] > 0:
                        while (samples[i+1] > 0) and i < (self.layer_count-1):
                            top = (x, y, start_z + i*SPACING, 20*samples[i])
                            bot = (x, y, start_z + (i+1)*SPACING, 20*samples[i+1])
                            self.cylinders.append((top, bot))
                            i += 1
                    elif samples[i] > 0:
                        self.num_circles += 1
                        raw_points.extend(self.create_circle(20*samples[i], x, y, start_z + i*SPACING))
                        colors.extend(leaf_color*360)
                    i += 1
                
            else:
                for i in range(self.layer_count):
                    self.num_circles += 1
                    raw_points.extend(self.create_circle(.001, x, y, start_z + i*SPACING))
                    colors.extend(leaf_color*360)
                    #colors.extend([0.4, 0.2, 0.2, 1.0]*360)
       

        #add the edge geometry to raw_points
        #TODO:again, incredibly inefficient 
        cur_z = start_z
        for i in range(self.layer_count):
            for e in self.edges:
                x1, y1, z1 = e.get_parent_coords()
                x2, y2, z2 = e.get_child_coords()
                raw_points.extend([x1, y1, cur_z, 1, x2, y2, cur_z, 1])
                colors.extend([0.4, 0.2, 0.2, 1.0]*2)
                #colors.extend([1.0, 0.0, 0.0, 1.0]*2)
                self.num_edges += 2
            cur_z += SPACING

        raw_points.extend(colors)
        #convert all geometry to a numpy array
        self.vertices = numpy.array(raw_points,
                                    dtype=numpy.float32)

        self.fragment_shader = None
        self.vertex_shader   = None
        self.shader_program  = None
        self.vao = None



    def key_check(self):
        '''
        Check rotation and zoom triggers, and
        act accordingly.
        '''
        if self.rot_y_left:
            self.y_deg += 1 
            glutPostRedisplay()
        elif self.rot_y_right:
            self.y_deg -= 1 
            glutPostRedisplay()
        if self.rot_x_up:
            self.x_deg += 1
            glutPostRedisplay()
        elif self.rot_x_down:
            self.x_deg -= 1
            glutPostRedisplay()
        if self.zoom_in:
            self.zoom_val -= 6 
            glutPostRedisplay()
        elif self.zoom_out:
            self.zoom_val += 6
            glutPostRedisplay()
        if self.rot_z_right:
            self.z_deg -= 1
            glutPostRedisplay()
        elif self.rot_z_left:
            self.z_deg += 1
            glutPostRedisplay()
            
    def special_key_press(self, key, x, y):
        '''
        Check for key press. If so, 
        set rotations to true. 
        '''
        if key == GLUT_KEY_LEFT:
            self.rot_y_left = 1
        elif key == GLUT_KEY_RIGHT:
            self.rot_y_right = 1
        elif key == GLUT_KEY_UP:
            self.rot_x_up = 1
        elif key == GLUT_KEY_DOWN:
            self.rot_x_down = 1
        glutPostRedisplay()

    def char_key_press(self, key, x, y):
        '''
        '''
        if key == 'c' or key == 'C':
            self.rot_z_right = 1
        elif key == 'z' or key == 'Z':
            self.rot_z_left = 1
        glutPostRedisplay()


    def special_key_release(self, key, x, y):
        '''
        Check for key release. This is used for 
        shutting down rotation signals. 
        '''
        if key == GLUT_KEY_LEFT:
            self.rot_y_left = 0
        elif key == GLUT_KEY_RIGHT:
            self.rot_y_right = 0
        elif key == GLUT_KEY_UP:
            self.rot_x_up = 0
        elif key == GLUT_KEY_DOWN:
            self.rot_x_down = 0
        glutPostRedisplay()

    def char_key_release(self, key, x, y):
        if key == 'c' or key == 'C':
            self.rot_z_right = 0
        elif key == 'z' or key == 'Z':
            self.rot_z_left = 0
        glutPostRedisplay()

        
    def mouse_button(self, button, state, x, y):
        '''
        Check for mouse clicks, and set the zoom trigger
        on and off accordingly. 
        '''
        if (button == GLUT_LEFT_BUTTON and state == GLUT_DOWN):
            self.zoom_in = 1
        elif (button == GLUT_LEFT_BUTTON and state == GLUT_UP):
            self.zoom_in = 0
        elif (button == GLUT_RIGHT_BUTTON and state == GLUT_DOWN):
            self.zoom_out = 1
        elif (button == GLUT_RIGHT_BUTTON and state == GLUT_UP):
            self.zoom_out = 0
    
    def display(self):
        '''
        Display the geometry. 
        '''

        #testing multi layer dipslay
        num_nodes = self.num_circles #len(self.nodes)*self.layer_count
        num_edges = len(self.edges)*self.layer_count

        glLoadIdentity()
        gluLookAt(0,0,self.leaf_count*15+self.start_zoom + self.zoom_val,
                  0,0,0,
                  0,1,0)
        glClearColor(1.0, 1.0, 1.0, 1.0)
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        glUseProgram(self.shader_program)
        glBindVertexArray(self.vao)

        glPushMatrix()
        glRotatef(self.x_deg, 1, 0, 0)
        glRotatef(self.y_deg, 0, 1, 0)
        glRotatef(self.z_deg, 0, 0, 1)
        glDisable(GL_CULL_FACE)

        #Draw the nodes
        i = 0
        while i < num_nodes:
            glDrawArrays(GL_TRIANGLE_FAN, i*360, 360)
            i += 1

        #Draw the branches
        start = i*360
        i = 0
        while i < (num_edges):
            glDrawArrays(GL_LINES, i*2 + start, 2)    
            i += 1

        #Draw cylinders -- this needs to be moved to the gpu
        quadric = gluNewQuadric()
        for c in self.cylinders:
            top = c[0]
            bot = c[1]
            top_x, top_y, top_z, top_r = top
            bot_x, bot_y, bot_z, bot_r = bot
            glPushMatrix()
            glTranslate(top_x, top_y, top_z) 
            gluCylinder(quadric, bot_r, top_r, abs(top_z-bot_z), 20, 20)
            glPopMatrix() 
       
        glPopMatrix()
        glBindVertexArray(0)
        glUseProgram(0)
        glutSwapBuffers()
      

    def create_circle(self, radius, x, y, z):
        deg2rad = 3.14159/180
        temp = [] 
        for i in range(0, 360):
            radian = i*deg2rad
            temp.append(x+math.cos(radian)*radius) 
            temp.append(y+math.sin(radian)*radius)
            temp.append(z)
            temp.append(1)

        return temp

    
    def execute(self):
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(1000,1000)
        glutCreateWindow('TreeViewer')

        #vec4( 0.4, 0.2, 0.2, 1 );
        #create the shaders
        self.fragment_shader = shaders.compileShader("""#version 130
        in  vec4 theColor;
        out vec4 fragColor;
        void main() {
            fragColor = theColor;
        }""", GL_FRAGMENT_SHADER)

        self.vertex_shader = shaders.compileShader("""#version 130
        in vec4 color;
        out vec4 theColor;
        void main() {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
            theColor    = color;
        }""", GL_VERTEX_SHADER)


        #previous way of creating the shader program
        #self.shader_program = shaders.compileProgram(self.vertex_shader,
        #                                             self.fragment_shader)


        self.shader_program = glCreateProgram()
        glBindAttribLocation(self.shader_program, 0, "position")
        glBindAttribLocation(self.shader_program, 1, "color")
        glAttachShader(self.shader_program, self.vertex_shader)
        glAttachShader(self.shader_program, self.fragment_shader)
        glLinkProgram(self.shader_program)

        #set up buffers on the gpu
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, vbo)
        glBufferData(GL_ARRAY_BUFFER, self.vertices.nbytes, self.vertices,
                     GL_STATIC_DRAW)
        glEnableVertexAttribArray(0)

        glEnableVertexAttribArray(1)
        #offset = 4*4*(self.num_circles*360 + len(self.edges)*2)
       
        offset = 4*4*(self.num_circles*360 + self.num_edges)

        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, None)
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, ctypes.c_void_p(offset))

        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0)

        #set up lighting and perspective
        glClearColor(1.,1.,1.,1.)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        lightZeroPosition = [10.,4.,10.,1.]
        lightZeroColor = [1.0,1.0,1.0,1.0] 
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
        glEnable(GL_LIGHT0)
        glutDisplayFunc(self.display)
        glMatrixMode(GL_PROJECTION)
        gluPerspective(1.,1.,-900.,900.)
        glMatrixMode(GL_MODELVIEW)

        #glutIgnoreKeyRepeat(1)
        glutSpecialFunc(self.special_key_press)
        glutKeyboardFunc(self.char_key_press)
        glutIgnoreKeyRepeat(1)
        glutSpecialUpFunc(self.special_key_release)
        glutKeyboardUpFunc(self.char_key_release)
        glutIdleFunc(self.key_check)
    
        glutMouseFunc(self.mouse_button)
        #glutMotionFunc(mouse_motion)

        glPushMatrix()
        glutMainLoop()


if __name__ == '__main__': 
    parser = argparse.ArgumentParser()
    parser.add_argument('newick_file', type=str)
    parser.add_argument('condensed_counts_file', type=str)
    parser.add_argument('layer_count', type=int, help="the number of samples to display",
                         nargs='?', default=MAX_LAYERS)
    args = parser.parse_args()
    newick_file = args.newick_file
    c_file      = args.condensed_counts_file
    layers      = (args.layer_count if args.layer_count <= MAX_LAYERS 
                  and args.layer_count > 0 else MAX_LAYERS)

    newick_f = open(newick_file, 'r')

    c_map    = CountsMap(c_file)
    newick_s = newick_f.readlines()[0]
    tree     = NewickTree(newick_s, c_map) 
    #if tree == 0:
    #   print('ERROR: failed to build tree')
    #   sys.exit()
    tv       = TreeViewer(tree.get_nodes(), tree.get_edges(), tree.total_leaves,
                          tree.get_num_samples(), layers) 
    tv.execute()

