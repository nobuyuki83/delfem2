import sys, cv2
sys.path.append("../module_py")
import dfm2

def main():
  path_img = "../test_inputs/lenna.png"
  np_img = cv2.imread(path_img)
  axis = dfm2.AxisXYZ(100)
  tex = dfm2.get_texture(np_img)
  dfm2.winDraw3d([axis,tex],winsize=(400,300),
                 camera_eye_up=[0,0,+1, 0,-1,0])

if __name__ == "__main__":
  main()