from OpenGL.GL import *
from OpenGL.GLUT import *
import dfm2


def print(x, y, font, text, color):
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()
  glColor3f(color[0], color[1], color[2])
  glRasterPos2f(x, y)
  for ch in text:
    glutBitmapCharacter(font, ctypes.c_int(ord(ch)))


def draw_sphere(pos, rad, color):
  if pos is None:
    return
  glColor3d(color[0], color[1], color[2])
  glTranslatef(+pos[0], +pos[1], +pos[2])
  glutSolidSphere(rad, 32, 32)
  glTranslatef(-pos[0], -pos[1], -pos[2])


class WindowManager:
  def __init__(self,view_height):
    self.camera = dfm2.gl.Camera(view_height)
    self.modifier = 0
    self.mouse_x = 0.0
    self.mouse_y = 0.0

  def special(self,key,x,y):
    if key == int(GLUT_KEY_PAGE_UP):
      self.camera.scale *= 1.03
    elif key == int(GLUT_KEY_PAGE_DOWN):
      self.camera.scale *= (1.0 / 1.03)
    glutPostRedisplay()

  def mouse(self,button, state, x, y):
    self.modifier = glutGetModifiers()
    self.mouse_x, self.mouse_y = dfm2.gl.mouse_screen_pos(x, y)
    glutPostRedisplay()

  def motion(self,x, y):  
    if self.modifier == 1: # shift
      self.mouse_x, self.mouse_y, self.camera.quat, self.camera.trans = dfm2.gl.motion_rot(
          x, y, self.mouse_x, self.mouse_y, self.camera.quat,
          self.camera.trans, self.camera.view_height)
    if self.modifier == 2 or self.modifier == 4: # ctrl or cmnd+alt
      self.mouse_x, self.mouse_y, self.camera.quat, self.camera.trans = dfm2.gl.motion_trans(
          x, y, self.mouse_x, self.mouse_y, self.camera.quat,
          self.camera.trans, self.camera.view_height)
    glutPostRedisplay()

class GlutWindow:
  def __init__(self,view_height):
    self.wm = WindowManager(view_height)
    self.draw_func = None

    glutInit()
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)  # zBuffer
    glutInitWindowSize(600, 600)
    glutInitWindowPosition(100, 100)
    glutCreateWindow("Visualizzatore_2.0")


  def keyboard(self, bkey, x, y):
    key = bkey.decode("utf-8")
    if key == 'q':
      exit()
    glutPostRedisplay()      

  def display(self):
    glClearColor(0.3, 0.5, 0.8, 1.0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glEnable(GL_DEPTH_TEST)  
    self.wm.camera.set_gl_camera()
    self.draw_func()
    glutSwapBuffers()  

  def reshape(self, width, height):
    glViewport(0, 0, width, height)

  def idle(self):
    glutPostRedisplay()

  def special(self,key, x, y):  
    self.wm.special(key,x,y)
    glutPostRedisplay()    

  def mouse(self,button, state, x, y):
    self.wm.mouse(button,state,x,y)
    glutPostRedisplay()    

  def motion(self,x, y):
    self.wm.motion(x,y)    
    glutPostRedisplay()  

  def draw_loop(self,draw_func0):  

    # Register callbacks
    glutReshapeFunc(self.reshape)
    glutDisplayFunc(self.display)
    glutMouseFunc(self.mouse)
    glutMotionFunc(self.motion)
    glutKeyboardFunc(self.keyboard)
    glutSpecialFunc(self.special)
    glutIdleFunc(self.idle)

    dfm2.dfm2.set_some_lighting()
    self.draw_func = draw_func0
    glutMainLoop()


