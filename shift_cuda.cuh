#ifndef _SHIFTE_KERNEL_
#define _SHIFTE_KERNEL_

#define WARP_SIZE 32
#define LOG_WARP_SIZE 5

//kernels for shift function
__global__ void 
count_leftright(int *rightcount, int *rightflag, int *leftcount, int *leftflag, 
                float *zpart, float pi, 
                float zeta1, float zeta0, float pi2_inv, int nparam, int mp, int m0)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int wid = tid >> LOG_WARP_SIZE;
  int lane = tid & (WARP_SIZE - 1);
	
  int bit;
  float zetaright, zetaleft;

  int right, left;
	
  //set the right and left flag
  if(tid<mp-m0+1)
  {
    zetaright=min(2.0*pi,zpart(3,tid+m0))-zeta1;
    zetaleft=zpart(3,tid+m0)-zeta0;
		
    if(zetaright*zetaleft>0)
    {
      zetaright=zetaright*pi2_inv;
      zetaright=zetaright-(float)(floorf(zetaright));
			
      if(zetaright<0.5) 
        right = 1, left = 0;
      else
        right = 0, left = 1;
    }
    else
    {
      right = 0, left = 0;
    }
    rightflag[tid] = right, leftflag[tid] = left;
  }

  //reduce the valid particles of each virtual processor
  if(tid<mp-m0+1)
  {
    bit = __ballot(right);
    if(lane==0)
      rightcount[wid]=__popc(bit);
    bit = __ballot(left);
    if(lane==0)
      leftcount[wid]=__popc(bit);
  }
	
}

__global__ void
compute_msendbuffer(int *msendbuffer, int *rightoffset, int *leftoffset, 
                    int *rightcount, int *leftcount, int mp, int m0)
{
  msendbuffer[0] = rightoffset[(mp-m0+1+WARPSIZE-1)/WARPSIZE-1] +
    rightcount[(mp-m0+1+WARPSIZE-1)/WARPSIZE-1];
  msendbuffer[1] = leftoffset[(mp-m0+1+WARPSIZE-1)/WARPSIZE-1] + 
    leftcount[(mp-m0+1+WARPSIZE-1)/WARPSIZE-1];
}

__global__ void
copy_right(float *sendright, int *hole_pos, 
           float *zpart, float *zpart0, 
           int *rightflag, const int *rightoffset, 
           int nzphase, int nparam, int mp, int mremain, int m0)
{
  int row_id, offset, bit, bit1;
	
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int wid = tid >> LOG_WARP_SIZE;
  int lane = tid & (WARP_SIZE - 1);
	
  if(tid<mp-m0+1)
  {
    offset = rightoffset[wid];         //the offset of each virtual processor
    int right = rightflag[tid];

    bit = __ballot(right);
		
    if (right == 1)
    {
      bit1 = bit & ((1<<lane)-1);
      for(row_id=1; row_id<=nparam; row_id++)
      {
        sendright(row_id, offset+(__popc(bit1))+1)=zpart(row_id, tid+m0);
        sendright(row_id+nparam, offset+(__popc(bit1))+1)=zpart0(row_id, tid+m0);
      }
      //remember the holes from particles going to right
      if(tid<mremain-m0+1)
        hole_pos[offset+(__popc(bit1))]=tid+m0; 
    }
  }
}

__global__ void
compute_righthole_count(int *righthole_count, int *rightflag, const int *rightoffset, int mremain, int m0)
{
  int bit;

  int wid = ((mremain - m0 + 1 + WARP_SIZE - 1) >> LOG_WARP_SIZE) - 1;
  int pid = wid*WARP_SIZE + threadIdx.x;

  if (pid < mremain - m0 + 1)
  {
    bit = __ballot(rightflag[pid]);
    if (threadIdx.x == 0)
      righthole_count[0] = __popc(bit) + rightoffset[wid];
  }   
}

__global__ void
copy_left(float *sendleft, int *hole_pos, 
          float *zpart, float *zpart0, 
          int *leftflag, int *leftoffset, 
          int nzphase, int nparam, int mp, int mremain, int m0)
{
  int row_id, offset, bit, bit1;
	
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int wid = tid >> LOG_WARP_SIZE;
  int lane = tid & (WARP_SIZE - 1);
	
  if(tid<mp-m0+1)
  {
    offset=leftoffset[wid];         //the offset of each virtual processor
    int left = leftflag[tid];
		
    bit = __ballot(left);

    if(left == 1)
    {
      bit1=bit & ((1<<lane)-1);
      for(row_id=1; row_id<=nparam; row_id++)
      {
        sendleft(row_id, offset+(__popc(bit1))+1)=zpart(row_id, tid+m0);
        sendleft(row_id+nparam, offset+(__popc(bit1))+1)=zpart0(row_id, tid+m0);
      }
      //remember the holes from particles going to left
      if(tid<mremain-m0+1)
        hole_pos[offset+(__popc(bit1))]=tid+m0;
    }
  }
}

__global__ void
copy_leftright(float *sendleft, float *sendright, 
               int *lefthole_pos, int *righthole_pos,
               float *zpart, float *zpart0, 
               int *leftflag, int *rightflag, 
               const int *leftoffset, const int *rightoffset, 
               int nzphase, int nparam, int mp, int mremain, int m0)
{
  int row_id, offset, bitl, bitr, bit1;
	
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int wid = tid >> LOG_WARP_SIZE;
  int lane = tid & (WARP_SIZE - 1);
	
  if(tid<mp-m0+1)
  {
    int left = leftflag[tid];
    int right = rightflag[tid];

    bitl = __ballot(left); 
    bitr = __ballot(right);

    if (left == 1)
    {
      offset = leftoffset[wid];         //the offset of each virtual processor
      bit1 = bitl & ((1<<lane)-1);
      for(row_id=1; row_id<=nparam; row_id++)
      {
        sendleft(row_id, offset+(__popc(bit1))+1)=zpart(row_id, tid+m0);
        sendleft(row_id+nparam, offset+(__popc(bit1))+1)=zpart0(row_id, tid+m0);
      }
      //remember the holes from particles going to left
      if(tid<mremain-m0+1)
        lefthole_pos[offset+(__popc(bit1))]=tid+m0; 
    }
		
    if (right == 1)
    {
      offset = rightoffset[wid];         //the offset of each virtual processor
      bit1 = bitr & ((1<<lane)-1);
      for(row_id=1; row_id<=nparam; row_id++)
      {
        sendright(row_id, offset+(__popc(bit1))+1)=zpart(row_id, tid+m0);
        sendright(row_id+nparam, offset+(__popc(bit1))+1)=zpart0(row_id, tid+m0);
      }
      //remember the holes from particles going to right
      if(tid<mremain-m0+1)
        righthole_pos[offset+(__popc(bit1))]=tid+m0; 
    }

  }
}

__global__ void
fill_hole(float *zpart, float *zpart0, 
           int *innerflag, int *hole_pos, int *innerpos, 
           int msend, int mremain, int m0, int nparam)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
		
  if(tid<msend)
  {
    int hpos = innerpos[tid];
    int inner = innerflag[tid];

    if (inner == 1) {
      int hole_id = hole_pos[hpos];
      for (int i = 1; i <= nparam; i++) {
        zpart(i,hole_id) = zpart(i,msend-1-tid + (mremain+1));
        zpart0(i,hole_id) = zpart0(i,msend-1-tid + (mremain+1));
      }
    }
  }
}


__global__ void
copy_recvleft(float *zpart, float *zpart0, 
              float *recvleft, int mrecvleft, 
              int nparam, int me, int nzphase)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;
	
  int m_row;
	
  if(tid<mrecvleft)
  {
    for(m_row=1; m_row<=nparam; m_row++)
    {
      zpart(m_row, tid+1+me)=recvleft(m_row, tid+1);
      zpart0(m_row, tid+1+me)=recvleft(m_row+nparam, tid+1);
    }
  }
}

__global__ void
copy_recvright(float *zpart, float *zpart0, 
               float *recvright, int mrecvleft, int mrecvright, 
               int nparam, int me, int nzphase)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;
	
  int m_row;
	
  if(tid<mrecvright)
  {
    for(m_row=1; m_row<=nparam; m_row++)
    {
      zpart(m_row, tid+1+me+mrecvleft)=recvright(m_row, tid+1);
      zpart0(m_row, tid+1+me+mrecvleft)=recvright(m_row+nparam, tid+1);
    }
  }
}

__global__ void compute_innerflag(int *innerflag, int *leftflag, int *rightflag, int n)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int tid2 = n - 1 - tid;
  if (tid < n)
  {
    innerflag[tid] = (leftflag[tid2]==0 && rightflag[tid2]==0) ? 1 : 0;
  }
}


__global__ void compute_sendoffset(int *sendoffset, int *leftoffset, int *rightoffset, int noffset)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if (tid < noffset) {
    sendoffset[tid] = leftoffset[tid] + rightoffset[tid];
  }
}

__global__ void compute_hole_pos(int *hole_pos, int *lefthole_pos, int *righthole_pos, int *offset, 
                                 int *leftflag, int *rightflag, int mp, int m0, int mremain)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int wid = tid >> LOG_WARP_SIZE;
  int lane = tid & (WARP_SIZE - 1);

  if (tid < mp-m0+1) {
    int warp_offset = offset[wid];
    unsigned int left = leftflag[tid];
    unsigned int right = rightflag[tid];
    
    int bit = __ballot(left | right);

    if (left == 1 || right == 1) {
      int bit1 = bit & ((1<<lane)-1);
      //if (tid < mremain-m0+1)
      hole_pos[warp_offset+(__popc(bit1))] = tid + m0;
    }
  }
}

#endif
