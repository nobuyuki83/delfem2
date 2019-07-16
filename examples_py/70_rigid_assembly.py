import sys
sys.path.append("..")
import PyDelFEM2 as dfm2
import PyDelFEM2.eigen as rb
import PyDelFEM2.gl._glfw

rb0 = rb.RigidBody(1.0, [+0.0, 1.0, +0.0])
rb1 = rb.RigidBody(0.1, [-1.0, 0.5, -1.0])
rb2 = rb.RigidBody(0.1, [-1.0, 0.5, +1.0])
rb3 = rb.RigidBody(0.1, [+1.0, 0.5, -1.0])
rb4 = rb.RigidBody(0.1, [+1.0, 0.5, +1.0])

rb1.add_contact_point([-1.1, 0, -1.1])
rb2.add_contact_point([-1.1, 0, +1.1])
rb3.add_contact_point([+1.1, 0, -1.1])
rb4.add_contact_point([+1.1, 0, +1.1])

jt0 = rb.Joint(0,1, [-1,+1,-1])
jt1 = rb.Joint(0,2, [-1,+1,+1])
jt2 = rb.Joint(0,3, [+1,+1,-1])
jt3 = rb.Joint(0,4, [+1,+1,+1])

rb_asm = rb.RigidBodyAssembly_Static([rb0,rb1,rb2,rb3,rb4],[jt0,jt1,jt2,jt3])
axis = dfm2.gl.AxisXYZ(1)
dfm2.gl._glfw.winDraw3d([axis,rb_asm],winsize=(400,300))
