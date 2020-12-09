#include <stdio.h>
#include <xc.h>

int main()
{
  xc_func_type func;

  double rho[5] = {0.11, 0.2, 0.3, 0.4, 0.5};
  double sigma[5] = {0.2, 0.3, 0.4, 0.5, 0.6};
  double exc[5];
  double vrho[5];
  double vsigma[5];
  const int N = 5;
  
  int i, vmajor, vminor, vmicro;
  int func_id = 130; // GGA_C_PBE

  xc_version(&vmajor, &vminor, &vmicro);
  printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

  if( xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0 ) {
    fprintf(stderr, "Functional '%d' not found\n", func_id);
    return 1;
  }

  xc_gga_exc_vxc(&func, 5, rho, sigma, exc, vrho, vsigma);

  for(i = 0; i < 5; i++) {
    printf("%18.10f %18.10f %18.10f %18.10f\n", rho[i], exc[i], vrho[i], vsigma[i]);
  }

  xc_func_end(&func);
}
