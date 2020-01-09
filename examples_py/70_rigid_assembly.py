
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

rb0 = dfm2.RigidBody(1.0, [+0.0, 1.0, +0.0])
rb1 = dfm2.RigidBody(0.1, [-1.0, 0.5, -1.0])
rb2 = dfm2.RigidBody(0.1, [-1.0, 0.5, +1.0])
rb3 = dfm2.RigidBody(0.1, [+1.0, 0.5, -1.0])
rb4 = dfm2.RigidBody(0.1, [+1.0, 0.5, +1.0])

rb1.add_contact_point([-1.1, 0, -1.1])
rb2.add_contact_point([-1.1, 0, +1.1])
rb3.add_contact_point([+1.1, 0, -1.1])
rb4.add_contact_point([+1.1, 0, +1.1])

jt0 = dfm2.Joint(0,1, [-1,+1,-1])
jt1 = dfm2.Joint(0,2, [-1,+1,+1])
jt2 = dfm2.Joint(0,3, [+1,+1,-1])
jt3 = dfm2.Joint(0,4, [+1,+1,+1])

rb_asm = dfm2.RigidBodyAssembly_Static([rb0,rb1,rb2,rb3,rb4],[jt0,jt1,jt2,jt3])
axis = dfm2.gl.AxisXYZ(1)
dfm2.gl.glfw.winDraw3d([axis,rb_asm],winsize=(400,300))
