#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include "cuda_util.h"
#define TF thrust::device_ptr<float>
#define TT thrust::device_ptr<int>


#define WARPSIZE  32
#define SCAN_BLOCK 1024

#define zpart(i, j)      	zpart[(j-1)*(nparam)+(i-1)]
#define zpart0(i, j)      	zpart0[(j-1)*(nparam)+(i-1)]
#define sendright(i, j)		sendright[(j-1)*nzphase+(i-1)]
#define sendleft(i, j)		sendleft[(j-1)*nzphase+(i-1)]
#define recvleft(i, j) 		recvleft[(j-1)*nzphase+(i-1)]
#define recvright(i, j)		recvright[(j-1)*nzphase+(i-1)]

#include "shift_cuda.cuh"

__global__ void myscan2(int *idata1, int *odata1, int *idata2, int *odata2, int n);
__global__ void myscan1(int *idata, int *odata, int n);

#define TIMING(x) 

#define MEM_FACTOR 2

extern "C" 
void shift_cuda(float *zpart, float *zpart0, int mpmax,int *mp,
    int nparam, int mtoroidal, int mype, int left_pe, int right_pe,
    int myrank_toroidal, int toroidal_comm0, int numberpe,
    float pi, float zeta0, float zeta1)
// mp needs to be pointer b/c it will be changed inside shifti
{
  if (numberpe == 1) return;

//   if (mype == 14) {
//     size_t freemem, totalmem;
//     cudaMemGetInfo(&freemem, &totalmem);
//     printf("P(%d): free = %f MB, total = %f MB\n",
//            mype, (double)freemem/pow(1024.0,2), (double)totalmem/pow(1024.0,2));
//   }

  MPI_Comm toroidal_comm = MPI_Comm_f2c(toroidal_comm0);
  float pi2_inv = 0.5/pi;
  TIMING( CTimer alltimer; CTimer timer );

  // static double tottime=0;
  // static double cpytime=0,cputime=0,mpitime=0,inittime=0,sendrecvtime=0;
  // static double keret[12]={0,0,0,0,0,0,0,0,0,0,0,0};
  // static size_t mpi_size=0;

  int mremain,m0,msend,msendtmp,mrecvtmp,idest,isource,isendtag,irecvtag,nzphase,
    isendcount,irecvcount,iteration;
  int msendleft,msendright,mrecvleft,mrecvright,msendbuffer[2];
	
  MPI_Status istatus;
  float *recvleft,*recvright,*sendleft,*sendright;
  TIMING( alltimer.Start() );

  TIMING( timer.Start() );
  static int p=(int)ceil((float)(mpmax)/WARPSIZE);
  static gpu_array<int> d_rightflag(mpmax);
  static gpu_array<int> d_rightcount(p);
  static gpu_array<int> d_rightoffset(p);
  static gpu_array<int> d_offset(p);
	
  static gpu_array<int> d_leftflag(mpmax);
  static gpu_array<int> d_leftcount(p);
  static gpu_array<int> d_leftoffset(p);
	
  static gpu_array<int> d_lefthole_pos(mpmax);
  static gpu_array<int> d_righthole_pos(mpmax);
  static gpu_array<int> d_hole_pos(mpmax);
  static gpu_array<int> d_msendbuffer(2);
  static gpu_array<int> d_righthole_count(1);
  static gpu_array<int> d_innerflag(mpmax);

  static float *d_sendright, *d_sendleft;
  static float *d_recvright, *d_recvleft;
  static int *d_innerpos;
  static int mem_sendright = 0,
    mem_sendleft = 0,
    mem_recvright = 0,
    mem_recvleft = 0,
    mem_send = 0;

  TIMING( inittime += timer.GetET(); timer.Start() ); 
  nzphase=2*(nparam);   //nzphase=14 if track_particles=1,=12 otherwise
  m0=1;
  iteration=0;
  TIMING( cputime += timer.GetET() );
shift: 
  //timer.Start();
  iteration=iteration+1;
	
  if(iteration>mtoroidal)
  {
    printf("endless particle sorting loop at PE=%d. \n",mype);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
	
  msend=0;
  msendright=0;
  msendleft=0;

  //Step 1. copy outgoing particles to send buffers
	
  //find the id of the outgoing electrons
  if(m0<=*mp)
  {	
    TIMING( timer.Start() );
    //find left and right particles using stream technology
    count_leftright<<<((*mp-m0+1)+127)/128,128>>>(d_rightcount.ptr(),d_rightflag.ptr(),
                                                  d_leftcount.ptr(),d_leftflag.ptr(),
                                                  zpart,pi,zeta1,zeta0,
                                                  pi2_inv,nparam,*mp,m0);
    TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[0] += timer.GetET(); timer.Start() );

    //compute the processor offsets
    if (((*mp-m0+1+WARPSIZE-1)/WARPSIZE) < 30000)
      myscan2<<<2,SCAN_BLOCK>>>(d_rightcount.ptr(),d_rightoffset.ptr(),
                                d_leftcount.ptr(),d_leftoffset.ptr(),
                                ((*mp-m0+1+WARPSIZE-1)/WARPSIZE));
    else 
    {
      thrust::exclusive_scan((TT)d_rightcount.ptr(),(TT)d_rightcount.ptr()+((*mp-m0+1+WARPSIZE-1)/WARPSIZE),
                             (TT)d_rightoffset.ptr());
      thrust::exclusive_scan((TT)d_leftcount.ptr(),(TT)d_leftcount.ptr()+((*mp-m0+1+WARPSIZE-1)/WARPSIZE),
                             (TT)d_leftoffset.ptr());
    }
    TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[1] += timer.GetET() ); 
    
    TIMING( timer.Start() );
    //get number of send particles  
    compute_msendbuffer<<<1,1>>>(d_msendbuffer.ptr(),d_rightoffset.ptr(),d_leftoffset.ptr(),
                                 d_rightcount.ptr(),d_leftcount.ptr(),*mp,m0);
    CUDA_SAFE_CALL( cudaGetLastError() );
    TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[2] += timer.GetET() ); 
    
    TIMING( timer.Start() );
    d_msendbuffer.D2H(msendbuffer,2);
    TIMING( cpytime += timer.GetET() );

    msendright=msendbuffer[0];
    msendleft =msendbuffer[1];
	
    //outgoing particles 
    msend=msendright+msendleft;
    // if (mype == 2) printf("pe=%d,send count:%d,%d,%d\n",mype,msend,msendright,msendleft);
  }

  if(msend!=(msendleft+msendright))
  {
    printf("mype=%d,msend NOT equal to msendleft+msendright",mype);
    msend=msendleft+msendright;	
  }

  if(iteration>1)
  {	
    //test: at least 1 particle needs to be shifted.
    mrecvtmp=0;
    msendtmp=0;
    if(msend>0)
      msendtmp=1;
    TIMING( timer.Start() );
    MPI_Allreduce(&msendtmp,&mrecvtmp,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    TIMING( mpitime += timer.GetET() ); 

    //no particle to be shifted,return
    if(mrecvtmp==0)
    {
      TIMING( tottime += alltimer.GetET() );
      return;
    }
  }
	
  //an extra space to prevent zero size when msendright(1)=msendleft(1)=0
  sendright=(float*)malloc(nzphase*max(1,msendright)*sizeof(float));
  sendleft=(float*)malloc(nzphase*max(1,msendleft)*sizeof(float));

  if (msendright >= mem_sendright) {
    if (mem_sendright != 0) 
      CUDA_SAFE_CALL(cudaFree(d_sendright));
    mem_sendright = MEM_FACTOR*msendright;
    CUDA_SAFE_CALL(cudaMalloc(&d_sendright, nzphase*max(1,mem_sendright)*sizeof(float)));
  } 
  if (msendleft >= mem_sendleft) {
    if (mem_sendleft != 0) 
      CUDA_SAFE_CALL(cudaFree(d_sendleft));
    mem_sendleft = MEM_FACTOR*msendleft;
    CUDA_SAFE_CALL(cudaMalloc(&d_sendleft,  nzphase*max(1,mem_sendleft)*sizeof(float)));
  }
  TIMING( inittime += timer.GetET() );

  //not outgoing particles
  mremain=*mp-msendright-msendleft;
  
  TIMING( timer.Start() );
  //copy right particles to sendright
  if (m0 <= *mp)
    copy_right<<<((*mp-m0+1)+255)/256,256>>>(d_sendright,d_righthole_pos.ptr(),
                                             zpart,zpart0,d_rightflag.ptr(),
                                             d_rightoffset.ptr(),nzphase,nparam,
                                             *mp,mremain,m0);
  CUDA_SAFE_CALL( cudaGetLastError() );
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[3] += timer.GetET() ); 
  TIMING( timer.Start() );
//   d_sendright.D2HAsync(sendright,d_sendright.size(),0);
  CUDA_SAFE_CALL( cudaMemcpy(sendright, d_sendright, nzphase*msendright*sizeof(float),
                             cudaMemcpyDeviceToHost) );
//   cudaEventRecord(event_d2h_right);
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); cpytime += timer.GetET() );

  //copy left particles to sendleft
  TIMING( timer.Start() );
  if (m0 <= *mp)
    copy_left<<<((*mp-m0+1)+255)/256,256>>>(d_sendleft,d_lefthole_pos.ptr(),
                                            zpart,zpart0,d_leftflag.ptr(),
                                            d_leftoffset.ptr(),nzphase,nparam, 
                                            *mp,mremain,m0);
  CUDA_SAFE_CALL( cudaGetLastError() );
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[5] += timer.GetET() );

  TIMING( timer.Start() );
//   d_sendleft.D2HAsync(sendleft,d_sendleft.size(),0);
  CUDA_SAFE_CALL( cudaMemcpy(sendleft, d_sendleft, nzphase*msendleft*sizeof(float),
                             cudaMemcpyDeviceToHost) );
//   cudaEventRecord(event_d2h_left);

  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() );  cpytime += timer.GetET() );
  // Step 2: fill in the hole of outgoing particle locations with non-outgoing particles	

//   gpu_array<int> d_innerpos(msend);
  if (msend >= mem_send) {
    if (mem_send != 0) 
      CUDA_SAFE_CALL(cudaFree(d_innerpos));
    mem_send = MEM_FACTOR*msend;
    CUDA_SAFE_CALL( cudaMalloc(&d_innerpos, max(1,mem_send)*sizeof(int)) );
  }

  TIMING( timer.Start() );
  if (m0 <= *mp) {
    int noffset = ((*mp-m0+1) + WARP_SIZE - 1) / WARP_SIZE;
    compute_sendoffset<<<(noffset+127)/128,128>>>
      (d_offset.ptr(), d_leftoffset.ptr(), d_rightoffset.ptr(), noffset);
    compute_hole_pos<<<((*mp-m0+1)+255)/256,256>>>
      (d_hole_pos.ptr(), d_lefthole_pos.ptr(), d_righthole_pos.ptr(), d_offset.ptr(),
       d_leftflag.ptr(), d_rightflag.ptr(), *mp, m0, mremain);
  }
	
  // find inner particles in *mremain+1:mp
  if (msend > 0)
    compute_innerflag<<<(msend+255)/256,256>>>(d_innerflag.ptr()+(mremain-m0+1),
                                               d_leftflag.ptr()+(mremain-m0+1),
                                               d_rightflag.ptr()+(mremain-m0+1),
                                               msend);
  CUDA_SAFE_CALL( cudaGetLastError() );
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[7] += timer.GetET() );
  
  TIMING( timer.Start() );
  // compute inner particle indices 
  if (msend < 10000)
    myscan1<<<1,SCAN_BLOCK>>>(d_innerflag.ptr() + (mremain-m0+1),d_innerpos,msend);
  else
    thrust::exclusive_scan((TT)d_innerflag.ptr()+(mremain-m0+1),
                           (TT)d_innerflag.ptr()+(*mp-m0+1),
                           (TT)d_innerpos);
  CUDA_SAFE_CALL( cudaGetLastError() );
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[8] += timer.GetET() );


  TIMING( timer.Start() );
  if (msend > 0)
    fill_hole<<<(msend+255)/256,256>>>(zpart,zpart0,
                                       d_innerflag.ptr()+(mremain-m0+1),
                                       d_hole_pos.ptr(),
                                       d_innerpos,msend,
                                       mremain,m0,nparam); 
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[8] += timer.GetET() );

  //send # of particle to move right
  mrecvleft=0;
  idest=right_pe;
  isource=left_pe;
  isendtag=myrank_toroidal;
  irecvtag=isource;

  TIMING( timer.Start() );
  MPI_Sendrecv(&msendright,1,MPI_INT,idest,isendtag,&mrecvleft,
               1,MPI_INT,isource,irecvtag,toroidal_comm,&istatus);
  TIMING( mpitime += timer.GetET() );
  TIMING( timer.Start() );
  recvleft=(float*)malloc(nzphase*max(1,mrecvleft)*sizeof(float));
  if (mrecvleft >= mem_recvleft) {
    if (mem_recvleft != 0) 
      CUDA_SAFE_CALL(cudaFree(d_recvleft));
    mem_recvleft = MEM_FACTOR*mrecvleft;
    CUDA_SAFE_CALL( cudaMalloc(&d_recvleft, nzphase*max(1,mem_recvleft)*sizeof(float)) );
  }
	
  //send particle to right and receive from left
  memset(recvleft,0,nzphase*max(1,mrecvleft)*sizeof(float));
  TIMING( inittime += timer.GetET() );
  isendcount=max(1,msendright)*nzphase;
  irecvcount=max(1,mrecvleft)*nzphase;
  // MPI wait for data transfer to complete
//   cudaEventSynchronize(event_d2h_right);
  TIMING( timer.Start() );
  MPI_Sendrecv(sendright,isendcount,MPI_FLOAT,idest,isendtag,recvleft,irecvcount,MPI_FLOAT,isource,irecvtag,toroidal_comm,&istatus);
  TIMING( sendrecvtime += timer.GetET() );
  TIMING( mpi_size += isendcount+irecvcount ); 
  TIMING( timer.Start() );
  cudaMemcpy(d_recvleft, recvleft, nzphase*mrecvleft*sizeof(float),
             cudaMemcpyHostToDevice);
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); cpytime += timer.GetET() );
  TIMING( timer.Start() );
  if (mrecvleft > 0)
    copy_recvleft<<<(mrecvleft+255)/256,256>>>(zpart,zpart0,
                                               d_recvleft,mrecvleft,
                                               nparam,mremain,nzphase);
  CUDA_SAFE_CALL( cudaGetLastError() );
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[10] += timer.GetET() );

  //send # of particle to move left
  mrecvright=0;
  idest=left_pe;
  isource=right_pe;
  isendtag=myrank_toroidal;
  irecvtag=isource;
  TIMING( timer.Start() );
  MPI_Sendrecv(&msendleft,1,MPI_INT,idest,isendtag,&mrecvright,1,MPI_INT,isource,irecvtag,toroidal_comm,&istatus);
  TIMING( mpitime += timer.GetET() );

  TIMING( timer.Start() );
  recvright=(float*)malloc(nzphase*max(1,mrecvright)*sizeof(float));
  if (mrecvright >= mem_recvright) {
    if (mem_recvright != 0)
      CUDA_SAFE_CALL(cudaFree(d_recvright));
    mem_recvright = MEM_FACTOR*mrecvright;
    CUDA_SAFE_CALL( cudaMalloc(&d_recvright, nzphase*max(1,mem_recvright)*sizeof(float)) );
  }
	
  //send particle to left and receive from right
  memset(recvright,0,nzphase*max(1,mrecvright)*sizeof(float));
  TIMING( inittime += timer.GetET() );  
  isendcount=max(1,msendleft)*nzphase;
  irecvcount=max(1,mrecvright)*nzphase;
  // MPI wait for data transfer to complete
//   cudaEventSynchronize(event_d2h_left);
  TIMING( timer.Start() );
  MPI_Sendrecv(sendleft,isendcount,MPI_FLOAT,idest,isendtag,recvright,irecvcount,MPI_FLOAT,isource,irecvtag,toroidal_comm,&istatus);
  TIMING( sendrecvtime += timer.GetET() );
  TIMING( mpi_size += isendcount+irecvcount ); 
  TIMING( timer.Start() );
//   d_recvright.H2DAsync(recvright,d_recvright.size(),0);
  cudaMemcpy(d_recvright, recvright, nzphase*mrecvright*sizeof(float),
             cudaMemcpyHostToDevice);
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); cpytime += timer.GetET() );
  TIMING( timer.Start() );
  if (mrecvright > 0) 
    copy_recvright<<<(mrecvright+255)/256,256>>>(zpart,zpart0,
                                                 d_recvright,mrecvleft,mrecvright,
                                                 nparam,mremain,nzphase);
  CUDA_SAFE_CALL( cudaGetLastError() );  
  TIMING( CUDA_SAFE_CALL( cudaDeviceSynchronize() ); keret[11] += timer.GetET() );

  //need extra particle array
  if(mremain+mrecvleft+mrecvright>mpmax)
  { 
    printf("need bigger particle array,%d,%d,%d. \n",mype,mpmax,mremain+mrecvleft+mrecvright);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
    
  //receive the electrons form others nodes
  *mp=mremain+mrecvleft+mrecvright;

  free(sendleft);
  free(sendright);
  free(recvleft);
  free(recvright);

  m0=*mp-mrecvright-mrecvleft+1;
	
  goto shift;
}

__forceinline__ __device__ uint scanwarp(uint val, volatile uint* sData, int maxlevel)
{
  // The following is the same as 2 * RadixSort::WARP_SIZE * warpId + threadInWarp =
  // 64*(threadIdx.x >> 5) + (threadIdx.x & (RadixSort::WARP_SIZE - 1))
  int localId = threadIdx.x;
  int idx = 2 * localId - (localId & (WARP_SIZE - 1));
  sData[idx] = 0;
  idx += WARP_SIZE;
  sData[idx] = val;

  if (0 <= maxlevel) { sData[idx] += sData[idx - 1]; }
  if (1 <= maxlevel) { sData[idx] += sData[idx - 2]; }
  if (2 <= maxlevel) { sData[idx] += sData[idx - 4]; }
  if (3 <= maxlevel) { sData[idx] += sData[idx - 8]; }
  if (4 <= maxlevel) { sData[idx] += sData[idx -16]; }

  return sData[idx] - val;  // convert inclusive -> exclusive
}

__forceinline__ __device__ uint4 scan_block(int *idata, volatile uint* ptr)
{

  uint idx = threadIdx.x;
  uint lane = idx & (WARP_SIZE - 1);
  uint wid = idx >> 5;

  uint4 val4 = *((uint4*)idata+idx);
  uint sum[3];
  sum[0] = val4.x;
  sum[1] = val4.y + sum[0];
  sum[2] = val4.z + sum[1];

  uint val = val4.w + sum[2];

  val = scanwarp(val, ptr, 4);
  __syncthreads();

  if (lane == WARP_SIZE - 1)
  {
    ptr[wid] = val + val4.w + sum[2];
  }
  __syncthreads();

  if (idx < WARP_SIZE)
    ptr[idx] = scanwarp(ptr[idx], ptr, 4);

  __syncthreads();

  val += ptr[wid];

  val4.x = val;
  val4.y = val + sum[0];
  val4.z = val + sum[1];
  val4.w = val + sum[2];

  return val4;

}
 
__global__ void myscan2(int *idata1, int *odata1, int *idata2, int *odata2, int n)
{
  int idx = threadIdx.x;
  __shared__ int input[4*SCAN_BLOCK];
  __shared__ uint ptr[2*SCAN_BLOCK];
  __shared__ uint tmp[1];

  int *idata, *odata;
  if (blockIdx.x == 0)
  {
    idata = idata1;
    odata = odata1;
  }
  if (blockIdx.x == 1)
  {
    idata = idata2;
    odata = odata2;
  }

  int niter = (int)ceil((float)n/(4*SCAN_BLOCK));

  for (int i = 0; i < niter; i++)
  {

    if (4*(idx+1) <= n)
    {
      for (int i = 0; i < 4; i++)
        input[4*threadIdx.x + i] = idata[4*idx+i]; 
    }
    else if (4*idx < n)
    {
      for (int i = 0; i < 4; i++)
        input[4*threadIdx.x+i] = (4*idx+i < n) ? idata[4*idx+i] : 0;
    }
    else
    {
      for (int i = 4*threadIdx.x; i < 4*(threadIdx.x+1); i++)
        input[i] = 0;
    }
    __syncthreads();

    uint4 val4;
    val4 = scan_block(input, ptr);

    if (i > 0)
    {
      val4.x += tmp[0];
      val4.y += tmp[0];
      val4.z += tmp[0];
      val4.w += tmp[0];
    }

    if (4*(idx+1) <= n)
    {
      odata[4*idx  ] = val4.x;
      odata[4*idx+1] = val4.y;
      odata[4*idx+2] = val4.z;
      odata[4*idx+3] = val4.w;
    }
    else if (4*idx < n)
    {
      odata[4*idx] = val4.x;
      if (4*idx + 1 < n)
        odata[4*idx+1] = val4.y;
      if (4*idx + 2 < n)
        odata[4*idx+2] = val4.z;
    }

    if (threadIdx.x == SCAN_BLOCK - 1 && i != niter - 1)
    {
      uint v = input[4*threadIdx.x + 3];
      tmp[0] = val4.w + v;
    }
    __syncthreads();

    idx += SCAN_BLOCK;
  }
}

__global__ void myscan1(int *idata, int *odata, int n)
{
  int idx = threadIdx.x;
  __shared__ int input[4*SCAN_BLOCK];
  __shared__ uint ptr[2*SCAN_BLOCK];
  __shared__ uint tmp[1];

  int niter = (int)ceil((float)n/(4*SCAN_BLOCK));

  for (int iter = 0; iter < niter; iter++)
  {
    if (4*(idx+1) <= n)
    {
      for (int i = 0; i < 4; i++)
        input[4*threadIdx.x + i] = idata[4*idx+i]; 
    }
    else if (4*idx < n)
    {
      for (int i = 0; i < 4; i++)
        input[4*threadIdx.x+i] = (4*idx+i < n) ? idata[4*idx+i] : 0;
    }
    else
    {
      for (int i = 4*threadIdx.x; i < 4*(threadIdx.x+1); i++)
        input[i] = 0;
    }
    __syncthreads();

    uint4 val4;
    val4 = scan_block(input, ptr);

    if (iter > 0)
    {
      val4.x += tmp[0];
      val4.y += tmp[0];
      val4.z += tmp[0];
      val4.w += tmp[0];
    }
    
    if (4*(idx+1) <= n)
    {
      odata[4*idx  ] = val4.x;
      odata[4*idx+1] = val4.y;
      odata[4*idx+2] = val4.z;
      odata[4*idx+3] = val4.w;
    }
    else if (4*idx < n)
    {
      odata[4*idx] = val4.x;
      if (4*idx + 1 < n)
        odata[4*idx+1] = val4.y;
      if (4*idx + 2 < n)
        odata[4*idx+2] = val4.z;
    }

    if (threadIdx.x == SCAN_BLOCK - 1 && iter != niter - 1)
    {
      uint v = input[4*threadIdx.x + 3];
      tmp[0] = val4.w + v;
    }
    __syncthreads();

    idx += SCAN_BLOCK;
  }
}

