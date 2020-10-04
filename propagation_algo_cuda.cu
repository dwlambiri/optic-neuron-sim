__global__  void propagation_algo_cuda(
        int N,
        float* tox_prod,
        short* axons,
        unsigned char* blue,
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
        float  outsideDetox,
        bool   algo,
        float* deathThr)
{
	//int idx = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
    int xidx = blockIdx.x * blockDim.x + threadIdx.x;
    int yidx = blockIdx.y * blockDim.y + threadIdx.y;

    int index = yidx*N+ xidx;
    int indexUp = (yidx-1)*N + xidx;
    int indexDown = (yidx+1)*N + xidx;
    int indexLeft = yidx*N + xidx -1;
    int indexRight = yidx*N + xidx+1;

    const short aliveAxon_c = 1;
    const short deadAxon_c  = 2;
    const short noAxon_c    = -1;

    const int stopProduction_c = 0;
    const unsigned char blueSky_c = 255;

    
    if((indexUp < lowerLimit) || (indexDown < lowerLimit) || (indexLeft < lowerLimit) || (indexRight < lowerLimit)) {
        return;
    }

    if((indexUp > upperLimit) || (indexDown > upperLimit) || (indexLeft > upperLimit) || (indexRight > upperLimit)) {
        return;
    }
   
    int centerIndex = centers[index];

    float extraAmount = 0;

    if(axons[index] == aliveAxon_c){
      if(cMap1[index] > deathThr[index]) {
        axons[index] = deadAxon_c;
        extraAmount = amountReleasedOnDeath;
        tox_prod[index] = stopProduction_c;
        detox[index] = outsideDetox;
        blue[index] = blueSky_c;
      }
    }
    
    if(centerIndex > 0 && axons[centerIndex-1] == deadAxon_c && tox_prod[index] > 0) {
        extraAmount = amountReleasedOnDeath;
        tox_prod[index] = stopProduction_c;
        detox[index] = outsideDetox;
        blue[index] = blueSky_c;
    }

    float di = dOutside;

    if(centerIndex > 0 && axons[centerIndex-1] == aliveAxon_c) {
       di = dInside; 
    }

    float t = cMap1[index];

    if(algo == true) {
        cMap2[index] = t +

                    (cMap1[indexUp] - t) * ((centers[indexUp]== -1)?0:di) +
                    (cMap1[indexDown] - t) * ((centers[indexDown]== -1)?0:di) +
                    (cMap1[indexLeft] - t) * ((centers[indexLeft]== -1)?0:di) +
                    (cMap1[indexRight] - t) * ((centers[indexRight]== -1)?0:di) +
                    tox_prod[index] + extraAmount;
    }
    else {
        cMap2[index] = t +
                (cMap1[indexUp] - t) * (di) +
                (cMap1[indexDown] - t) * (di) +
                (cMap1[indexLeft] - t) * (di) +
                (cMap1[indexRight] - t) * (di) +
                 tox_prod[index] + extraAmount;  
    }

    cMap2[index] *= detox[index];
               
/*
    if(cMap2[index] > deathThreashold && tox_prod[index] > 0) {
        cMap2[index] = amountReleasedOnDeath;
        tox_prod[index] = 0; 
        return;
    }
*/      
    
}
