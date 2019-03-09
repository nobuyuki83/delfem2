#ifndef tex_h
#define tex_h

class CTexture
{
public:
  std::vector<unsigned char> aRGB;
  unsigned int id_tex;
  int h,w;
  
public:
  CTexture(){
    id_tex = 0;
    this->w = 0;
    this->h = 0;
    aRGB.clear();
  }
  CTexture(int w, int h, const unsigned char* pD, const std::string& typeData){
    this->h = h;
    this->w = w;
    this->aRGB.assign(pD,pD+h*w*3);
    if( typeData == "bgr" ){
      for(int i=0;i<h*w;++i){ // rgb -> bgr
        unsigned char b0 = aRGB[i*3+0];
        unsigned char r0 = aRGB[i*3+2];
        aRGB[i*3+0] = r0;
        aRGB[i*3+2] = b0;
      }
    }
    id_tex = 0;
  }
  
  void LoadTex()
  {
    if( id_tex == 0 ){
      ::glGenTextures(1, &id_tex);
    }
    glBindTexture(GL_TEXTURE_2D, id_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    assert( (int)aRGB.size() == w*h*3 );
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                 w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
                 aRGB.data() );
    glBindTexture(GL_TEXTURE_2D, 0);
  }
  
  void Draw(){
    if( id_tex == 0 ){ return; }
    /*
     const CVector3& dx = x_axis;
     const CVector3& dy = Cross(z_axis,dx);
     const double lx = lengrid*nResX;
     const double ly = lengrid*nResY;
     CVector3 p0 = origin;
     CVector3 p1 = origin + lx*dx;
     CVector3 p2 = origin + lx*dx + ly*dy;
     CVector3 p3 = origin + ly*dy;
     */
    ::glEnable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glBindTexture(GL_TEXTURE_2D, id_tex);
    ::glColor3d(1,1,1);
    ::glBegin(GL_QUADS);
    ::glTexCoord2d(0.0, 0.0); ::glVertex3d(0,0,0);
    ::glTexCoord2d(1.0, 0.0); ::glVertex3d(w,0,0);
    ::glTexCoord2d(1.0, 1.0); ::glVertex3d(w,h,0);
    ::glTexCoord2d(0.0, 1.0); ::glVertex3d(0,h,0);
    ::glEnd();
    ::glBindTexture(GL_TEXTURE_2D, 0);
    ::glDisable(GL_TEXTURE_2D);
  }
  
  std::vector<double>  MinMaxXYZ(){
    std::vector<double> m(6,0.0);
    m[0] = 0;
    m[1] = w;
    m[2] = 0;
    m[3] = h;
    m[4] = 0;
    m[5] = 0;
    return m;
  }
};

#endif
