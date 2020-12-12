#include <stdio.h>
#include <xc.h>

int main()
{
  xc_func_type func;

  double rho[1] = {1.1};
  double sigma[1] = {2.0};
  double lapl[1] = {0.0};
  double tau[1] = {1.0};

  // outputs
  double exc[1];
  double vrho[1];
  double vsigma[1];
  double vlapl[1];
  double vtau[1];
  
  const int N = 1;
  
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
    printf("%18.10f %18.10f: %18.10f %18.10f %18.10f %18.10f %18.10f\n",
        rho[i], sigma[i], 
        exc[i], vrho[i], vsigma[i], vlapl[i], vtau[i]);
  }

  xc_func_end(&func);
}
