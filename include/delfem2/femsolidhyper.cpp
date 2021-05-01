

#include "delfem2/femutil.h"
#include "delfem2/femsolidhyper.h"
#include <cmath>

void MakeStressConstitute_Hyper2D(
    const double c1,
    const double c2,
    const double dudx[][2],
    const double press,
    double stress_2ndpk[][2],
    double constit[][2][2][2],
    double &rcg_pinv1,
    double &rcg_pinv2,
    double &rcg_pinv3,
    double rcg_inv[][2])
{
  constexpr unsigned int ndim = 2;
  double strain_gl[ndim][ndim];
  for (unsigned int idim = 0; idim < ndim; idim++) {
    for (unsigned int jdim = 0; jdim < ndim; jdim++) {
      strain_gl[idim][jdim] = 0.5 * (dudx[idim][jdim] + dudx[jdim][idim]);
      for (unsigned int kdim = 0; kdim < ndim; kdim++) {
        strain_gl[idim][jdim] += 0.5 * dudx[kdim][idim] * dudx[kdim][jdim];
      }
    }
  }

  double rcg[ndim][ndim];
  for (unsigned int idim = 0; idim < ndim; idim++) {
    for (unsigned int jdim = 0; jdim < ndim; jdim++) {
      rcg[idim][jdim] = dudx[idim][jdim] + dudx[jdim][idim];
      for (unsigned int kdim = 0; kdim < ndim; kdim++) {
        rcg[idim][jdim] += dudx[kdim][idim] * dudx[kdim][jdim];
      }
    }
    rcg[idim][idim] += 1.0;
  }

  rcg_pinv1 = rcg[0][0] + rcg[1][1] + 1.0;
  rcg_pinv2 = rcg[0][0] * rcg[1][1] + rcg[0][0] + rcg[1][1] - rcg[0][1] * rcg[1][0];
  rcg_pinv3 = rcg[0][0] * rcg[1][1] - rcg[0][1] * rcg[1][0];
  {
    const double inv_pinv3 = 1.0 / rcg_pinv3;
    rcg_inv[0][0] = +inv_pinv3 * rcg[1][1];
    rcg_inv[0][1] = -inv_pinv3 * rcg[0][1];
    rcg_inv[1][0] = -inv_pinv3 * rcg[1][0];
    rcg_inv[1][1] = +inv_pinv3 * rcg[0][0];
  }

  const double inv13_rcg_pinv3 = 1.0 / pow(rcg_pinv3, 1.0 / 3.0);
  const double inv23_rcg_pinv3 = 1.0 / pow(rcg_pinv3, 2.0 / 3.0);

  {

    for (unsigned int i = 0; i < ndim * ndim; i++) { *(&stress_2ndpk[0][0] + i) = 0.0; }
    for (unsigned int idim = 0; idim < ndim; idim++) {
      for (unsigned int jdim = 0; jdim < ndim; jdim++) {
        stress_2ndpk[idim][jdim] -=
            2.0 * c2 * inv23_rcg_pinv3 * rcg[idim][jdim]
            + 2.0 * (c1 * rcg_pinv1 * inv13_rcg_pinv3 + c2 * 2.0 * rcg_pinv2 * inv23_rcg_pinv3) / 3.0
              * rcg_inv[idim][jdim];
      }
    }
    {
      double dtmp1 = 2.0 * c1 * inv13_rcg_pinv3 + 2.0 * c2 * inv23_rcg_pinv3 * rcg_pinv1;
      for (unsigned int idim = 0; idim < ndim; idim++) {
        stress_2ndpk[idim][idim] += dtmp1;
      }
    }
    {
      double dtmp1 = 2.0 * press * rcg_pinv3;
      for (unsigned int idim = 0; idim < ndim; idim++) {
        for (unsigned int jdim = 0; jdim < ndim; jdim++) {
          stress_2ndpk[idim][jdim] += dtmp1 * rcg_inv[idim][jdim];
        }
      }
    }
  }

  {
    for (unsigned int i = 0; i < ndim * ndim * ndim * ndim; i++) { *(&constit[0][0][0][0] + i) = 0.0; }
    for (unsigned int idim = 0; idim < ndim; idim++) {
      for (unsigned int jdim = 0; jdim < ndim; jdim++) {
        for (unsigned int kdim = 0; kdim < ndim; kdim++) {
          for (unsigned int ldim = 0; ldim < ndim; ldim++) {
            constit[idim][jdim][kdim][ldim]
                += 4.0 * c1 * inv13_rcg_pinv3 / 3.0 * (
                rcg_inv[idim][jdim] * rcg_inv[kdim][ldim] * rcg_pinv1 / 3.0
                + rcg_inv[idim][kdim] * rcg_inv[ldim][jdim] * rcg_pinv1 * 0.5
                + rcg_inv[idim][ldim] * rcg_inv[kdim][jdim] * rcg_pinv1 * 0.5);
            constit[idim][jdim][kdim][ldim]
                += 4.0 * c2 * inv23_rcg_pinv3 * 2.0 / 3.0 * (
                rcg_inv[idim][jdim] * rcg_inv[kdim][ldim] * rcg_pinv2 * (2.0 / 3.0)
                + rcg_inv[idim][jdim] * rcg[kdim][ldim]
                + rcg[idim][jdim] * rcg_inv[kdim][ldim]
                + rcg_inv[idim][kdim] * rcg_inv[jdim][ldim] * rcg_pinv2 * 0.5
                + rcg_inv[idim][ldim] * rcg_inv[jdim][kdim] * rcg_pinv2 * 0.5);
          }
        }
        double dtmp1 = 4.0 * c1 * inv13_rcg_pinv3 / 3.0 * rcg_inv[idim][jdim];
        for (unsigned int kdim = 0; kdim < ndim; kdim++) {
          constit[idim][jdim][kdim][kdim] -= dtmp1;
          constit[kdim][kdim][idim][jdim] -= dtmp1;
        }
        double dtmp2 = 4.0 * c2 * inv23_rcg_pinv3 * rcg_pinv1 * (2.0 / 3.0) * rcg_inv[idim][jdim];
        for (unsigned int kdim = 0; kdim < ndim; kdim++) {
          constit[idim][jdim][kdim][kdim] -= dtmp2;
          constit[kdim][kdim][idim][jdim] -= dtmp2;
        }
        constit[idim][idim][jdim][jdim] += 4.0 * c2 * inv23_rcg_pinv3;
        constit[idim][jdim][jdim][idim] -= 2.0 * c2 * inv23_rcg_pinv3;
        constit[idim][jdim][idim][jdim] -= 2.0 * c2 * inv23_rcg_pinv3;
      }
    }
    for (unsigned int idim = 0; idim < ndim; idim++) {
      for (unsigned int jdim = 0; jdim < ndim; jdim++) {
        for (unsigned int kdim = 0; kdim < ndim; kdim++) {
          for (unsigned int ldim = 0; ldim < ndim; ldim++) {
            constit[idim][jdim][kdim][ldim] +=
                4.0 * press * rcg_pinv3 *
                rcg_inv[idim][jdim] * rcg_inv[kdim][ldim]
                - 2.0 * press * rcg_pinv3 * (
                    rcg_inv[idim][kdim] * rcg_inv[ldim][jdim]
                    + rcg_inv[idim][ldim] * rcg_inv[kdim][jdim]);
          }
        }
      }
    }
  }
}

namespace delfem2 {
namespace femsolidhyper {

double WdWdCddWddC_Solid3Compression(
    double dWdC2[],
    double ddWddC2[][6],
    //
    const double dudx[][3])
{
  constexpr unsigned int ndim = 3;

  // right cauchy-green tensor
  double C[ndim][ndim];
  RightCauchyGreen_DispGrad<3>(C,dudx);

  double Cinv[ndim][ndim];
  const double p3C = DetInv_Mat3(Cinv,C);

  { // extracting independent component of symmetric tensor
    const double tmp0 = (p3C-1) * 2 * p3C;
    dWdC2[0] = tmp0 * Cinv[0][0];
    dWdC2[1] = tmp0 * Cinv[1][1];
    dWdC2[2] = tmp0 * Cinv[2][2];
    dWdC2[3] = tmp0 * Cinv[0][1];
    dWdC2[4] = tmp0 * Cinv[1][2];
    dWdC2[5] = tmp0 * Cinv[2][0];
  }
  { // Extracting independent components in the constitutive tensor
    constexpr unsigned int istdim2ij[6][2] = {
        {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {2, 0}
    };
    for (unsigned int istdim = 0; istdim < 6; istdim++) {
      for (unsigned int jstdim = 0; jstdim < 6; jstdim++) {
        const unsigned int idim = istdim2ij[istdim][0];
        const unsigned int jdim = istdim2ij[istdim][1];
        const unsigned int kdim = istdim2ij[jstdim][0];
        const unsigned int ldim = istdim2ij[jstdim][1];
        const double ddp3CddC = p3C * (
            +4.*Cinv[idim][jdim]*Cinv[kdim][ldim]
            -2.*Cinv[idim][ldim]*Cinv[jdim][kdim]
            -2.*Cinv[idim][kdim]*Cinv[jdim][ldim] );
        const double dp3CdpC_dp3CdpC = 4. * p3C * p3C * Cinv[idim][jdim] * Cinv[kdim][ldim];
        ddWddC2[istdim][jstdim] = (p3C-1.)*ddp3CddC + dp3CdpC_dp3CdpC;
      }
    }
  }
  return 0.5*(p3C-1)*(p3C-1);
}

/**
 * compute energy density and its gradient & hessian w.r.t. right Cauchy-Green tensor given the displacement gradient tensor
 * @param[out] dWdC2 gradient of energy density w.r.t. right Cauchy-Green tensor
 * @param[out] ddWddC2 hessian of energy density w.r.t. right Cauchy-Green tensor
 * @param[in] c1
 * @param[in] c2
 * @param[in] dudx displacement gradient tensor
 * @return density of elastic potential energy
 */
double WdWdCddWddC_Solid3HyperMooneyrivlin2Reduced(
    double dWdC2[],
    double ddWddC2[][6],
    //
    const double c1,
    const double c2,
    const double dudx[][3])
{
  constexpr unsigned int ndim = 3;
  
  // right cauchy-green tensor
  double C[ndim][ndim];
  RightCauchyGreen_DispGrad<3>(C,dudx);

  // invariants of Cauchy-Green tensor
  const double p1C = C[0][0] + C[1][1] + C[2][2];
  const double p2C =
      + C[0][0] * C[1][1]
      + C[0][0] * C[2][2]
      + C[1][1] * C[2][2]
      - C[0][1] * C[1][0]
      - C[0][2] * C[2][0]
      - C[1][2] * C[2][1];

  double Cinv[ndim][ndim];
  const double p3C = DetInv_Mat3(Cinv,C);

  const double tmp1 = 1.0 / pow(p3C, 1.0 / 3.0);
  const double tmp2 = 1.0 / pow(p3C, 2.0 / 3.0);
  const double pi1C = p1C*tmp1; // 1st reduced invariant
  const double pi2C = p2C*tmp2; // 2nd reduced invariant
  const double W = c1*(pi1C-3.) + c2*(pi2C-3.);

  { // compute 2nd Piola-Kirchhoff tensor here
    double S[ndim][ndim]; // 2nd Piola-Kirchhoff tensor
    for (unsigned int idim = 0; idim < ndim; idim++) {
      for (unsigned int jdim = 0; jdim < ndim; jdim++) {
        S[idim][jdim] =
            - 2.0 * c2 * tmp2 * C[idim][jdim]
            - 2.0 * (c1 * pi1C + c2 * 2.0 * pi2C) / 3.0 * Cinv[idim][jdim];
      }
    }
    {
      const double dtmp1 = 2.0 * c1 * tmp1 + 2.0 * c2 * tmp2 * p1C;
      S[0][0] += dtmp1;
      S[1][1] += dtmp1;
      S[2][2] += dtmp1;
    }
    { // 2nd piola-kirchhoff tensor is symmetric. Here extracting 6 independent elements.
      dWdC2[0] = S[0][0];
      dWdC2[1] = S[1][1];
      dWdC2[2] = S[2][2];
      dWdC2[3] = S[0][1];
      dWdC2[4] = S[1][2];
      dWdC2[5] = S[2][0];
    }
  }
  
  // computing constituive tensor from here
  double ddWddC[ndim][ndim][ndim][ndim];
  for (unsigned int i = 0; i < ndim * ndim * ndim * ndim; i++) { *(&ddWddC[0][0][0][0] + i) = 0.0; }
  for (unsigned int idim = 0; idim < ndim; idim++) {
    for (unsigned int jdim = 0; jdim < ndim; jdim++) {
      for (unsigned int kdim = 0; kdim < ndim; kdim++) {
        for (unsigned int ldim = 0; ldim < ndim; ldim++) {
          double tmp = 0;
          tmp += 4.0 * c1 * tmp1 / 3.0 * (
              Cinv[idim][jdim] * Cinv[kdim][ldim] * p1C / 3.0
              + Cinv[idim][kdim] * Cinv[ldim][jdim] * p1C * 0.5
              + Cinv[idim][ldim] * Cinv[kdim][jdim] * p1C * 0.5);
          tmp += 4.0 * c2 * tmp2 * 2.0 / 3.0 * (
              Cinv[idim][jdim] * Cinv[kdim][ldim] * p2C * (2.0 / 3.0)
              + Cinv[idim][jdim] * C[kdim][ldim]
              + C[idim][jdim] * Cinv[kdim][ldim]
              + Cinv[idim][kdim] * Cinv[jdim][ldim] * p2C * 0.5
              + Cinv[idim][ldim] * Cinv[jdim][kdim] * p2C * 0.5);
          ddWddC[idim][jdim][kdim][ldim] += tmp;
        }
      }
    }
  }
  for (unsigned int idim = 0; idim < ndim; idim++) {
    for (unsigned int jdim = 0; jdim < ndim; jdim++) {
      double dtmp1 = 4.0 * c1 * tmp1 / 3.0 * Cinv[idim][jdim];
      for (unsigned int kdim = 0; kdim < ndim; kdim++) {
        ddWddC[idim][jdim][kdim][kdim] -= dtmp1;
        ddWddC[kdim][kdim][idim][jdim] -= dtmp1;
      }
      const double dtmp2 = 4.0 * c2 * tmp2 * p1C * (2.0 / 3.0) * Cinv[idim][jdim];
      for (unsigned int kdim = 0; kdim < ndim; kdim++) {
        ddWddC[idim][jdim][kdim][kdim] -= dtmp2;
        ddWddC[kdim][kdim][idim][jdim] -= dtmp2;
      }
      ddWddC[idim][idim][jdim][jdim] += 4.0 * c2 * tmp2;
      ddWddC[idim][jdim][jdim][idim] -= 2.0 * c2 * tmp2;
      ddWddC[idim][jdim][idim][jdim] -= 2.0 * c2 * tmp2;
    }
  }
  { // Extracting independent components in the constitutive tensor
    constexpr unsigned int istdim2ij[6][2] = {
        {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {2, 0}
    };
    for (unsigned int istdim = 0; istdim < 6; istdim++) {
      for (unsigned int jstdim = 0; jstdim < 6; jstdim++) {
        const unsigned int idim = istdim2ij[istdim][0];
        const unsigned int jdim = istdim2ij[istdim][1];
        const unsigned int kdim = istdim2ij[jstdim][0];
        const unsigned int ldim = istdim2ij[jstdim][1];
        ddWddC2[istdim][jstdim] = ddWddC[idim][jdim][kdim][ldim];
      }
    }
  }
  return W;
}

void AdddWddW_EnergyGradHessian(
    double dW[8][3],
    double ddW[8][8][3][3],
    const double dudx[3][3],
    const double dndx[8][3],
    const double dWdC2[6],
    const double ddWddC2[6][6],
    double detwei)
{
  double dC2dU[8][3][6];
  {
    double z_mat[3][3];
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){
        z_mat[idim][jdim] = dudx[idim][jdim];
      }
      z_mat[idim][idim] += 1.0;
    }
    for(unsigned int ino=0;ino<8;ino++){
      for(unsigned int idim=0;idim<3;idim++){
        dC2dU[ino][idim][0] = dndx[ino][0]*z_mat[idim][0];
        dC2dU[ino][idim][1] = dndx[ino][1]*z_mat[idim][1];
        dC2dU[ino][idim][2] = dndx[ino][2]*z_mat[idim][2];
        dC2dU[ino][idim][3] = dndx[ino][0]*z_mat[idim][1]+dndx[ino][1]*z_mat[idim][0];
        dC2dU[ino][idim][4] = dndx[ino][1]*z_mat[idim][2]+dndx[ino][2]*z_mat[idim][1];
        dC2dU[ino][idim][5] = dndx[ino][2]*z_mat[idim][0]+dndx[ino][0]*z_mat[idim][2];
      }
    }
  }
  // make ddW
  for(unsigned int ino=0;ino<8;ino++){
    for(unsigned int jno=0;jno<8;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          double dtmp1 = 0.0;
          for(unsigned int gstdim=0;gstdim<6;gstdim++){
            for(unsigned int hstdim=0;hstdim<6;hstdim++){
              dtmp1 += ddWddC2[gstdim][hstdim]
                       *dC2dU[ino][idim][gstdim]*dC2dU[jno][jdim][hstdim];
            }
          }
          ddW[ino][jno][idim][jdim] += detwei*dtmp1;
        }
      }
      {
        double dtmp2 = 0.0;
        dtmp2 += dWdC2[0]*dndx[ino][0]*dndx[jno][0];
        dtmp2 += dWdC2[1]*dndx[ino][1]*dndx[jno][1];
        dtmp2 += dWdC2[2]*dndx[ino][2]*dndx[jno][2];
        dtmp2 += dWdC2[3]*(dndx[ino][0]*dndx[jno][1]+dndx[ino][1]*dndx[jno][0]);
        dtmp2 += dWdC2[4]*(dndx[ino][1]*dndx[jno][2]+dndx[ino][2]*dndx[jno][1]);
        dtmp2 += dWdC2[5]*(dndx[ino][2]*dndx[jno][0]+dndx[ino][0]*dndx[jno][2]);
        for(unsigned int idim=0;idim<3;idim++){
          ddW[ino][jno][idim][idim] += detwei*dtmp2;
        }
      }
    }
  }
  // make dW
  for(unsigned int ino=0;ino<8;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      double dtmp1 = 0.0;
      for(unsigned int istdim=0;istdim<6;istdim++){
        dtmp1 += dC2dU[ino][idim][istdim]*dWdC2[istdim];
      }
      dW[ino][idim] += detwei*dtmp1;
    }
  }
}

}
}



void delfem2::AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
    double& W,
    double dW[8][3],
    double ddW[8][8][3][3],
    double& vol,
    double c1,
    double c2,
    const double aP0[8][3],
    const double aU[8][3],
    unsigned int iGauss)
{
  vol = 0.0;
  const unsigned int nInt = delfem2::NIntLineGauss[iGauss];
  for(unsigned int ir1=0;ir1<nInt;ir1++) {
    for (unsigned int ir2 = 0; ir2 < nInt; ir2++) {
      for (unsigned int ir3 = 0; ir3 < nInt; ir3++) {
        double dndx[8][3];
        const double detwei = DiffShapeFuncAtQuadraturePoint_Hex(
            dndx,
            iGauss, ir1, ir2, ir3, aP0);
        vol += detwei;
        double dudx[3][3]; DispGrad_GradshapeDisp<3,8>( dudx, dndx, aU );
        double dWdC2[6], ddWddC2[6][6];
        const double w0 = femsolidhyper::WdWdCddWddC_Solid3HyperMooneyrivlin2Reduced(
            dWdC2, ddWddC2,
            c1, c2, dudx);
        W += w0*detwei;
        femsolidhyper::AdddWddW_EnergyGradHessian(dW, ddW,
            dudx, dndx, dWdC2, ddWddC2, detwei);
      } // r1
    } // r2
  } // r3
}


void delfem2::AddWdWddW_Solid3Compression_Hex(
    double& W,
    double dW[8][3],
    double ddW[8][8][3][3],
    double& vol,
    double stiff_comp,
    const double aP0[8][3],
    const double aU[8][3],
    unsigned int iGauss)
{
  vol = 0.0;
  const unsigned int nInt = delfem2::NIntLineGauss[iGauss];
  for(unsigned int ir1=0;ir1<nInt;ir1++) {
    for (unsigned int ir2 = 0; ir2 < nInt; ir2++) {
      for (unsigned int ir3 = 0; ir3 < nInt; ir3++) {
        double dndx[8][3];
        const double detwei = DiffShapeFuncAtQuadraturePoint_Hex(
            dndx,
            iGauss, ir1, ir2, ir3, aP0);
        vol += detwei;
        double dudx[3][3]; DispGrad_GradshapeDisp<3,8>( dudx, dndx, aU );
        double dWdC2[6], ddWddC2[6][6];
        const double w0 = femsolidhyper::WdWdCddWddC_Solid3Compression(
            dWdC2, ddWddC2,
            dudx);
        W += w0*detwei*stiff_comp;
        femsolidhyper::AdddWddW_EnergyGradHessian(
            dW, ddW,
            dudx, dndx, dWdC2, ddWddC2, detwei*stiff_comp);
      } // r1
    } // r2
  } // r3
}