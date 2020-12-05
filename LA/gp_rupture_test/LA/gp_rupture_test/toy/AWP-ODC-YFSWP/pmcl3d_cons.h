#define BLOCK_SIZE_X 2
#define BLOCK_SIZE_Y 2
#define BLOCK_SIZE_Z 4
#define align 32
#define loop  1 
#define ngsl 4     /* number of ghost cells x loop */
#define ngsl2 8  /* ngsl * 2 */

#define Both  0
#define Left  1
#define Right 2
#define Front 3
#define Back  4

#define NEDZ_EP 160 /*max k to save final plastic strain*/

#define CUCHK(call) {                                    \
  cudaError_t err = call;                                                    \
  if( cudaSuccess != err) {                                                \
  fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
          __FILE__, __LINE__, cudaGetErrorString( err) );              \
  fflush(stderr); \
  exit(EXIT_FAILURE);                                                  \
  } }

//precompiles variables for DM
#define MAXGRIDS 10

/*order in which variables are stored in swap buffers for transition zone
DO NOT CHANGE */
#define sbvpos_u1 0
#define sbvpos_v1 1
#define sbvpos_w1 2
#define sbvpos_xx 3
#define sbvpos_yy 4
#define sbvpos_zz 5
#define sbvpos_xy 6
#define sbvpos_xz 7
#define sbvpos_yz 8

//WEDMI window length = 2*WWL + 1
#define WWL 2

//HIGHEST ORDER OF FILTER is MAXFILT-1
#define MAXFILT 20

#define MPIRANKIO 400000
