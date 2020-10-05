__global__  void propagation_algo_cuda3d(
        int N,
        float* tox_prod,
        short* axons,
        unsigned char* blue,
        float* cMapResult,
        const float* cMapUp,
        const float* cMap,
        const float* cMapDown,
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
        float* deathThr,
        bool top,
        bool bottom,
        bool injury)
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
      if(top == true && cMap[index] > deathThr[index]) {
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

    float t = cMap[index];

    if(algo == true) {
        cMapResult[index] = t +
                    (cMap[indexUp] - t) * ((centers[indexUp]== -1)?0:di) +
                    (cMap[indexDown] - t) * ((centers[indexDown]== -1)?0:di) +
                    (cMap[indexLeft] - t) * ((centers[indexLeft]== -1)?0:di) +
                    (cMap[indexRight] - t) * ((centers[indexRight]== -1)?0:di);

                    if((top == false) && (bottom == false)) {
                        cMapResult[index] += (cMapUp[index] - t)*di +
                                             (cMapDown[index] -t)*di;
                    }
                    if((top == true) && (bottom == false)) {
                        cMapResult[index] += (cMapDown[index] -t)*di;
                    }
                    if((top == false) && (bottom == true)) {
                        cMapResult[index] += (cMapUp[index] -t)*di;
                    }
                    if(injury == true) {
                        cMapResult[index] += tox_prod[index] + extraAmount;
                    }
    }
    else {
        cMapResult[index] = t +
                (cMap[indexUp] - t) * (di) +
                (cMap[indexDown] - t) * (di) +
                (cMap[indexLeft] - t) * (di) +
                (cMap[indexRight] - t) * (di); 
        if((top == false) && (bottom == false)) {
            cMapResult[index] += (cMapUp[index] - t)*di +
                    (cMapDown[index] -t)*di;
        }
        if((top == true) && (bottom == false)) {
            cMapResult[index] += (cMapDown[index] -t)*di;
        }
        if((top == false) && (bottom == true)) {
            cMapResult[index] += (cMapUp[index] -t)*di;
        }
        if(injury == true) {
            cMapResult[index] += tox_prod[index] + extraAmount;
        }
    }

    cMapResult[index] *= detox[index];
               
/*
    if(cMapResult[index] > deathThreashold && tox_prod[index] > 0) {
        cMapResult[index] = amountReleasedOnDeath;
        tox_prod[index] = 0; 
        return;
    }
*/      
    
}
