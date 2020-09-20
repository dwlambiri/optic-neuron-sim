__global__  void propagation_algo_cuda(
        int N,
        float* tox_prod,
        short* axons,
        float* blue,
        float* cMap2,
        const float* cMap1,
        float* detox,
        const float* centers,
        float dInside,
        float dOutside,
        int lowerLimit,
        int upperLimit,
        float deathThreashold,
        float amountReleasedOnDeath,
        float  outsideDetox)
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
   
    int centerIndex = centers[index];

    if(axons[index] == 1){
      if(cMap2[index] > deathThreashold) {
        cMap2[index] = amountReleasedOnDeath;
        tox_prod[index] = 0; 
        axons[index] = 2;
        blue[index] = 64;
        detox[index] = outsideDetox;
        return;
      }
    }

    if(centerIndex > 0 && axons[centerIndex-1] ==2 && tox_prod[index] > 0) {
        cMap2[index] = amountReleasedOnDeath;
        tox_prod[index] = 0;
        detox[index] = outsideDetox;
        blue[index] = 64;
        return;
    }

    float di = dOutside;

    if(centerIndex > 0 && axons[centerIndex-1] == 1) {
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
    if(cMap2[index] > deathThreashold && tox_prod[index] > 0) {
        cMap2[index] = amountReleasedOnDeath;
        tox_prod[index] = 0; 
        return;
    }
*/      
    
}
