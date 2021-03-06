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
from pylab import *
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
    def __init__(self, tree, layers):
    
        self.layer_count = layers
        self.num_samples = tree.get_num_samples()
        self.num_leaves  = tree.get_num_leaves()
        self.nodes       = tree.get_nodes()
        self.edges       = tree.get_edges()
        self.radius      = tree.get_radius()
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
        self.plate_count = 0
        self.rim_count   = 0

        #Create the circle geometry
        start_z = -1*SPACING*(self.layer_count/2) #(self.num_samples/100)*-5
        raw_points = []
        self.num_circles = 0
        self.num_edge_points = 0

        colors = []
        cyl_colors = []
        cylinders  = [] 
        self.num_cylinders = 0

        color_map  = matplotlib.cm.get_cmap('rainbow')        
        leaf_count = 0
        
        #Prototype for leaf interpolation. 
        #TODO: we could probably make this signficantly more effecient
        for node in self.nodes:
            #FIXME: the leaf counter is still off for layered display
            x, y, z = node.get_coords()
            if node.is_leaf():
                samples = node.get_counts_list() 
                i = 0
                 
                #if we only have one layer, just create spheres. 
                if self.layer_count == 1:
                    self.num_circles += 1
                    raw_points.extend(self.create_circle(20*samples[i], x, y, start_z + i*SPACING))
                    colors.extend(color_map(float(leaf_count)/float(self.num_leaves))*360)

                #if we have multiple layers, look for cylinders.  
                while i < (self.layer_count-1):
                    if samples[i] > 0 and samples[i+1] > 0:
                        while (samples[i+1] > 0) and i < (self.layer_count-1):
                            cylinders.extend(self.create_cylinder(20*samples[i], 20*samples[i+1],
                                              x, y, start_z + i*SPACING, SPACING))
                            cyl_colors.extend(color_map(float(leaf_count)/float(self.num_leaves))*722)
                            self.num_cylinders += 1
                            i += 1
                    elif samples[i] > 0:
                        self.num_circles += 1
                        raw_points.extend(self.create_circle(20*samples[i], x, y, start_z + i*SPACING))
                        colors.extend(color_map(float(leaf_count)/float(self.num_leaves))*360)
                    i += 1
                leaf_count += 1
                
            else:
                for i in range(self.layer_count):
                    self.num_circles += 1
                    raw_points.extend(self.create_circle(.001, x, y, start_z + i*SPACING))
                    colors.extend([0.4, 0.2, 0.2, 1.0]*360)
       
        #Previous method of creating a full tree for every layer
        '''
        cur_z = start_z
        for i in range(self.layer_count):
            for e in self.edges:
                x1, y1, z1 = e.get_parent_coords()
                x2, y2, z2 = e.get_child_coords()
                raw_points.extend([x1, y1, cur_z+.2, 1, x2, y2, cur_z+.2, 1])
                raw_points.extend([x1, y1, cur_z-.2, 1, x2, y2, cur_z-.2, 1])
                #colors.extend([0.4, 0.2, 0.2, 1.0]*4)
                colors.extend([1.0, 0.9, 0.41, 1.0]*4)
                self.num_edge_points += 4
            cur_z += SPACING
        '''

        #create the tree branches
        cur_z = start_z
        if self.layer_count > 1:
            for i in range(2):
                for e in self.edges:
                    x1, y1, z1 = e.get_parent_coords()
                    x2, y2, z2 = e.get_child_coords()
                    raw_points.extend([x1, y1, cur_z+.2, 1, x2, y2, cur_z+.2, 1])
                    raw_points.extend([x1, y1, cur_z-.2, 1, x2, y2, cur_z-.2, 1])
                    #colors.extend([0.4, 0.2, 0.2, 1.0]*4)
                    colors.extend([1.0, 0.9, 0.41, 1.0]*4)
                    self.num_edge_points += 4
                cur_z = (start_z + (self.layer_count - 1)*SPACING)
       
        else:
            for e in self.edges:
                x1, y1, z1 = e.get_parent_coords()
                x2, y2, z2 = e.get_child_coords()
                raw_points.extend([x1, y1, cur_z+.2, 1, x2, y2, cur_z+.2, 1])
                raw_points.extend([x1, y1, cur_z-.2, 1, x2, y2, cur_z-.2, 1])
                #colors.extend([0.4, 0.2, 0.2, 1.0]*4)
                colors.extend([1.0, 0.9, 0.41, 1.0]*4)
                self.num_edge_points += 4


        raw_points.extend(cylinders)

        plate_colors = []
        plate_rims   = []
        rim_colors   = []
        
        #Create the plates TODO: this can probably go in the
        #branch loop
        cur_z = start_z
        if self.layer_count > 1:
            for i in range(2): 
                raw_points.extend(self.create_circle(self.radius-.1, 0, 0, cur_z))
                plate_colors.extend([0.86, 0.92, 0.95, 1.0]*360)
                self.plate_count += 1
                cur_z = (start_z + (self.layer_count - 1)*SPACING)

        else:
            raw_points.extend(self.create_circle(self.radius-.1, 0, 0, cur_z))
            plate_colors.extend([0.86, 0.92, 0.95, 1.0]*360)
            self.plate_count += 1
            
        #Create the layer rings/rims
        for i in range(self.layer_count):
            plate_rims.extend(self.create_circle(self.radius, 0, 0, start_z + i*SPACING))
            rim_colors.extend([0.0, 0.0, 0.0, 1.0]*360)
            self.rim_count += 1
           
        #old method for creating plates and rims for every layer
        '''
        for i in range(self.layer_count):
            raw_points.extend(self.create_circle(self.radius-.1, 0, 0, start_z + i*SPACING))
            plate_rims.extend(self.create_circle(self.radius, 0, 0, start_z + i*SPACING))
            #plate_colors.extend([0.77, 0.77, 0.77, 1.0]*360)
            plate_colors.extend([0.86, 0.92, 0.95, 1.0]*360)
            rim_colors.extend([0.0, 0.0, 0.0, 1.0]*360)
        '''

        #add the plate and rim geometry and colors to raw_points
        raw_points.extend(plate_rims)
        raw_points.extend(colors)
        raw_points.extend(cyl_colors)
        raw_points.extend(plate_colors)
        raw_points.extend(rim_colors)

        #convert all geometry to a numpy array
        self.vertices = numpy.array(raw_points,
                                    dtype=numpy.float32)

        self.fragment_shader = None
        self.vertex_shader   = None
        self.shader_program  = None
        self.vao             = None



    def key_check(self):
        '''
        Check rotation and zoom triggers, and
        act accordingly.
        '''
        if self.rot_y_left:
            self.y_deg -= 1 
            glutPostRedisplay()
        elif self.rot_y_right:
            self.y_deg += 1 
            glutPostRedisplay()
        if self.rot_x_up:
            self.x_deg += 1
            glutPostRedisplay()
        elif self.rot_x_down:
            self.x_deg -= 1
            glutPostRedisplay()
        if self.zoom_in:
            self.zoom_val += 1 
            glutPostRedisplay()
        elif self.zoom_out:
            self.zoom_val -= 1
            glutPostRedisplay()
        if self.rot_z_right:
            self.z_deg += 1
            glutPostRedisplay()
        elif self.rot_z_left:
            self.z_deg -= 1
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
        Check for c or z key presses. These keys
        are corresponding to rotation around the 
        z axis.  
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
        '''
        Check for c or z key release. These keys
        are corresponding to rotation around the 
        z axis.  
        '''
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

    def reshape(w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0)
        glMatrixMade(GL_MODELVIEW)
    
    def display(self):
        '''
        Display the geometry. 
        '''
        glLoadIdentity()
        gluLookAt(0,0, -1.0*self.num_leaves + self.start_zoom + self.zoom_val,
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


        num_nodes = self.num_circles #len(self.nodes)*self.layer_count
        #num_edges = len(self.edges)*self.layer_count
        num_edges = self.num_edge_points/2

        i = 0
        while i < num_nodes:
            glDrawArrays(GL_TRIANGLE_FAN, i*360, 360)
            i += 1

        #Draw the branches
        edge_start = i*360
        i = 0
        while i < num_edges:
            glDrawArrays(GL_LINES, i*2 + edge_start, 2)    
            i += 1

        #Draw cylinders
        cyl_start = i*2 + edge_start
        i = 0
        while i < self.num_cylinders:
            glDrawArrays(GL_TRIANGLE_STRIP, i*722 + cyl_start, 722)
            i += 1

        #Draw plates
        plate_start = i*722 + cyl_start
        i = 0
        while i < self.plate_count:
            glDrawArrays(GL_TRIANGLE_FAN, i*360 + plate_start, 360)
            i += 1

        #Draw plate rims
        rim_start = plate_start + i*360
        i = 0
        while i < self.rim_count:
            glDrawArrays(GL_LINE_LOOP, i*360 + rim_start, 360)
            i += 1

        glPopMatrix()
        glBindVertexArray(0)
        glUseProgram(0)
        glutSwapBuffers()
      

    def create_circle(self, radius, x, y, z):
        '''
        Create the geometry for a single circle. 
        An array filled with the geometry is returned. 
        Each point consists of x, y, z, 1
        '''
        deg2rad = math.pi/180.0
        points = [] 
        for i in range(0, 360):
            radian = float(i)*deg2rad
            points.append(x + math.cos(radian)*radius) 
            points.append(y + math.sin(radian)*radius)
            points.append(z)
            points.append(1)
        return points


    def create_cylinder(self, top_radius, bot_radius, x, y, z, h):
        '''
        Create the geometry for a cylinder and store this 
        geometry in an array that is returned. Each point
        is of the form x, y, z, 1
        '''
        deg2rad = math.pi/180.0
        points = []
        for i in range(361):
            radian = float(i)*deg2rad
            points.append(x + math.cos(radian)*top_radius)
            points.append(y + math.sin(radian)*top_radius)
            points.append(z)
            points.append(1.0)
            points.append(x + math.cos(radian)*bot_radius)
            points.append(y + math.sin(radian)*bot_radius)
            points.append(z + h)
            points.append(1.0)
        return points             

    #TODO: this should probably be called 'set-up' or something
    #      along those lines.  
    def execute(self):
        glutInit(sys.argv)

        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(1000,1000)
        glutCreateWindow('TreeViewer')
        win_w = glutGet(GLUT_WINDOW_WIDTH)
        win_h = glutGet(GLUT_WINDOW_HEIGHT)

        #create the shaders
        self.fragment_shader = shaders.compileShader("""#version 130
        in  vec4 theColor;
        out vec4 fragColor;
        uniform float go;
        void main() {
            if (go == 1.0)
                fragColor = vec4(1.0, 0.0, 0.0, 1.0);
            else
                fragColor = theColor;
        }""", GL_FRAGMENT_SHADER)

        self.vertex_shader = shaders.compileShader("""#version 130
        in vec4 color;
        out vec4 theColor;
        void main() {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
            theColor    = color;
        }""", GL_VERTEX_SHADER)

        #create the shader program
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
       
        color_offset = 4*4*(self.num_circles*360 + self.num_edge_points + self.num_cylinders*722 +
                            (self.plate_count+self.rim_count)*360)

        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, None)
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, ctypes.c_void_p(color_offset))

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
        #glutReshapeFunc(reshape)#TODO: for some reason, reshaping seems 
                                 #      to work but throws a value error.
                                 #      I'm not using it for now.

        w = glutGet(GLUT_WINDOW_WIDTH)
        h = glutGet(GLUT_WINDOW_HEIGHT)
        glViewport(0, 0, w, h)
        glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, 1000000000.0)
        glMatrixMode(GL_PROJECTION)
        gluPerspective(25.,float(w)/float(h),1.,1000000000.)
        glMatrixMode(GL_MODELVIEW)

        glutSpecialFunc(self.special_key_press)
        glutKeyboardFunc(self.char_key_press)
        glutIgnoreKeyRepeat(1)
        glutSpecialUpFunc(self.special_key_release)
        glutKeyboardUpFunc(self.char_key_release)
        glutIdleFunc(self.key_check)
    
        glutMouseFunc(self.mouse_button)

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
    tv       = TreeViewer(tree, layers) 
    tv.execute()

