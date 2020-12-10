#include <stdio.h>
#include <xc.h>

int main()
{
  xc_func_type func;

  double rho[5] = {1.01, 1.02, 0.03, 1.04, 1.05};
  double sigma[5] = {1.01, 0.02, 0.03, 0.04, 0.05};
  double lapl[5] = {2.1, 2.2, 2.3, 2.4, 2.6};
  double tau[5] = {1.1, 2.1, 3.1, 4.1, 5.1};

  // outputs
  double exc[5];
  double vrho[5];
  double vsigma[5];
  double vlapl[5];
  double vtau[5];
  
  const int N = 5;
  
  int i, vmajor, vminor, vmicro;
  int func_id = 263; // MGGA_X_SCAN
  //int func_id = 267; // MGGA_C_SCAN

  xc_version(&vmajor, &vminor, &vmicro);
  printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

  if( xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0 ) {
    fprintf(stderr, "Functional '%d' not found\n", func_id);
    return 1;
  }

  xc_mgga_exc_vxc(&func, N, rho, sigma, lapl, tau,
    exc, vrho, vsigma, vlapl, vtau);

  for(i = 0; i < N; i++) {
    printf("%18.10f: %18.10f %18.10f %18.10f %18.10f %18.10f\n",
        rho[i],
        exc[i], vrho[i], vsigma[i], vlapl[i], vtau[i]);
  }

  xc_func_end(&func);
}
