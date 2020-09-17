__global__  void propagation_algo_cuda(
        int N,
        float* tox_prod,
        short* axons,
        float* blue,
        float* cMap2,
        const float* cMap1,
        const float* detox,
        const float* centers,
        double dInside,
        double dOutside,
        int lowerLimit,
        int upperLimit)
{
	//int idx = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
    int xidx = blockIdx.x * blockDim.x + threadIdx.x;
    int yidx = blockIdx.y * blockDim.y + threadIdx.y;

    int index = yidx*N+ xidx;
    int indexUp = (yidx-1)*N + xidx;
    int indexDown = (yidx+1)*N + xidx;
    int indexLeft = yidx*N + xidx -1;
    int indexRight = yidx*N + xidx+1;

    
    if((indexUp < lowerLimit) || (indexDown < lowerLimit) || (indexLeft < lowerLimit) || (indexRight < lowerLimit)) {
        return;
    }

    if((indexUp > upperLimit) || (indexDown > upperLimit) || (indexLeft > upperLimit) || (indexRight > upperLimit)) {
        return;
    }
   
    if(centers[index] < 0) {
         return;
    }

    float di = dOutside;

    if(centers[index] > 0 ) {
       di = dInside; 
    }

    float t = cMap1[index];
    cMap2[index] = t +
                (cMap1[indexUp] - t) * di +
                (cMap1[indexDown] - t) * di +
                (cMap1[indexLeft] - t) * di +
                (cMap1[indexRight] - t) * di +
                tox_prod[index];

    cMap2[index] *= detox[index];
/*
    if(cMap2[index] > 22 && tox_prod[index] > 0) {
        cMap2[index] = 10000;
        tox_prod[index] = 0; 
        return;
    }
*/
       
    if(axons[index] == 1){
      if(cMap2[index] > 22) {
        cMap2[index] = 10000;
        tox_prod[index] = 0; 
        axons[index] = 2;
        blue[index] = 1;
      }
      return;
    }
    
    int centerIndex = centers[index];

    if(centerIndex <= 0){
        return;
    }
       
    if(tox_prod[index] > 0 && axons[centerIndex-1] == 2 && cMap2[index] > 22) {
        cMap2[index] = 10000;
        tox_prod[index] = 0;
        blue[index] = 1;
    }
}
