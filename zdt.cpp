#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
int iseed;
float tauV[800];
float tauE[800][800];
int nbVertices;
int nbEdges;
int degree[800];
std::vector<int> succ[800];
std::vector<int> nonSucc[800];
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
float randFloat() {
    int il,ih,it;
    double rc;
    ih = iseed/iq;
    il = iseed%iq;
    it = ia*il-ir*ih;
    if (it > 0) iseed = it;
    else iseed = ic+it;
    rc = ic;
    return iseed/rc;
}

int getNextRand(int limit) {
    return (int)(100000*randFloat()) % limit;
}

void init(){
    std::memset(tauV,0,sizeof(tauV));
    std::memset(tauE,0,sizeof(tauE));
    std::memset(degree,0,sizeof(degree));
    for(int i=0;i<800;++i){
        succ[i].clear();
        nonSucc[i].clear();
    }
}
void createGraph(){
    int i, j, nbSucc, nbNonSucc;
    int cp, src, dest;
    int edgesCounter=0;
    for(i=0;i<nbEdges;++i){
        scanf("%d%d",&src,&dest);
        src--,dest--;
        edgesCounter++;
        tauE[src][dest]=1.0;
        tauE[dest][src]=1.0;
        degree[src]++;
        degree[dest]++;
    }
    for(src=0;src<nbVertices;++src){
        dest=0;
        for(dest=0;dest<nbVertices;++dest){
            if(tauE[src][dest]>0){
                succ[src].push_back(dest);
            }else if(src!=dest){
                nonSucc[src].push_back(dest);
            }
        }
    }
}
void initPhero(float tauMax){
    int i,j;
    for(i = 0;i<nbVertices;++i){
        tauV[i]=tauMax;
    }
}
int choose(float *p,int nbCand){
    float f = (float) randFloat();
    int left=0;
    int right=nbCand-1;
    int k=0;
    float total = p[nbCand-1];
    while(left<right){
        k = (left+right+1)/2;
        if(f<(p[k-1]/total))
            right=k-1;
        else if(f>(p[k]/total))
            left = k+1;
        else
            return k;
    }
    if(left>=0 && left<nbCand)
        return left;
    return 0;
}
int isEdge(int i,int j){
    return tauE[i][j]>0;
}
int walk(int alpha,int *clique,int firstVertex){
    int i, v;
    int candidates[nbVertices];
    int nbCandidates;
    float total;
    float p[nbVertices];
    int cliqueSize = 1;
    
    clique[0]=firstVertex;
    nbCandidates = degree[firstVertex];
    total=0;
    for(i=0;i<succ[firstVertex].size();++i){
        candidates[i] = succ[firstVertex][i];
        p[i] = std::pow(tauV[candidates[i]],alpha) + total;
        total = p[i];
    }
    while(nbCandidates>0){
        i=choose(p,nbCandidates);
        v = candidates[i];
        clique[cliqueSize] = v;
        cliqueSize++;
        nbCandidates--;
        candidates[i] = candidates[nbCandidates];
        i = 0;
        total=0;
        while(i<nbCandidates){
            if(isEdge(v,candidates[i])){
                p[i] = std::pow(tauV[candidates[i]],alpha)+total;
                total = p[i];
                i++;
            }else{
                nbCandidates--;
                candidates[i]=candidates[nbCandidates];
            }
        }
    }
    return cliqueSize;
}

int repair(int* clique,int &cliqueSize){
    int vi, vj, vk, i, j, k, l, stop, found;
    int delta = 0;
    i = getNextRand(cliqueSize);
    stop = i;
    while(1){
        vi = clique[i];
        found = 0;
        for(j=0;((j<nonSucc[vi].size())&&(found==0));++j){
            vj = nonSucc[vi][j];
            for(k=j+1;((k<nonSucc[vi].size())&&(found==0));++k){
                vk = nonSucc[vi][k];
                if(isEdge(vj,vk)){
                    for(l=0;((l<cliqueSize)&&
                                ((clique[l]==vi)||
                                 ((isEdge(vj,clique[l]))&&(isEdge(vk,clique[l])))));
                            l++);
                    if(l==cliqueSize)
                        found=1;
                    else if(!isEdge(vj,clique[l]))
                        k=nbVertices;
                }
            }
        }
        if(found==1){
            clique[i]=vj;
            clique[cliqueSize]=vk;
            cliqueSize++;
            delta++;
            stop=i;
        }else{
            ++i;
            if(i==cliqueSize)
                i=0;
            if(i==stop)
                return delta;
        }
    }
}
int displayClique(int cliqueSize, int* clique){
    int i;
    sort(clique,clique+cliqueSize);
    printf("%d\n",cliqueSize);
    //printf("\t- clique size: %d vertices\n\t- clique list: ",cliqueSize);
    for(i=0 ; i<cliqueSize ; i++) 
        printf("%d ", clique[i]+1);
    printf("\n");
    return 1;
}

void evaporate(float rho,float tauMin){
    int i,j;
    for(i=0;i<nbVertices;++i){
        tauV[i] *= rho;
        if(tauV[i]<tauMin)
            tauV[i]=tauMin;
    }
}

void reinforce(int* clique,int cliqueSize,float qty,float tauMax){
    int i,j;
    for(i=0;i<cliqueSize;++i){
        tauV[clique[i]]+=qty;
        if(tauV[clique[i]]>tauMax)
            tauV[clique[i]]=tauMax;
    }
}

int main(){
    int alpha = 1;         
    int best = 10000;      
    float rho = 0.01;      
    int maxCycles = 5000;  
    int nbAnts = 30;       
    float tauMin = 0.01;   
    float tauMax = 6;      
    int mustRepair = 1;    
    iseed = 3;

    while(scanf("%d%d",&nbVertices,&nbEdges)==2){
        init();
        createGraph();
        initPhero(tauMax);
        if(best>nbVertices)
            best = nbVertices;
        int bestCliqueSize = 0;
        int bestClique[nbVertices];
        int currentCliqueSize;
        int currentClique[nbVertices];
        int bestCliqueCycleSize;
        int bestCliqueCycle[nbVertices];
        int nbRepair = 0;
        float qty;
        int i, j;
        int nbIter = 0;
        int nbCycles;
        for(nbCycles=0;((nbCycles<maxCycles)&&(bestCliqueSize<best));++nbCycles){
            bestCliqueCycleSize=0;
            for(i=0;((i<nbAnts)&&(bestCliqueSize<best));++i){
                currentCliqueSize = walk(alpha,currentClique,getNextRand(nbVertices));
                nbIter+=currentCliqueSize;
                if (currentCliqueSize > bestCliqueCycleSize) {
                    for (j=0; j<currentCliqueSize; j++) 
                        bestCliqueCycle[j] = currentClique[j];
                    bestCliqueCycleSize = currentCliqueSize;
                }
            }
            if((mustRepair==1)&&(bestCliqueCycleSize<best))
                nbRepair+= repair(bestCliqueCycle,bestCliqueCycleSize);
            if (bestCliqueCycleSize > bestCliqueSize) {
                bestCliqueSize=bestCliqueCycleSize;
                for (j=0; j<bestCliqueSize; j++) 
                    bestClique[j]=bestCliqueCycle[j];
                //printf("Best clique size = %d at cycle %d\n", bestCliqueSize, nbCycles);
            }
            evaporate(1-rho,tauMin);
            qty = 1.0/(float)(bestCliqueSize+1-bestCliqueCycleSize);
            reinforce(bestCliqueCycle,bestCliqueCycleSize,qty,tauMax);
        }
        displayClique(bestCliqueSize, bestClique);
    }

    return 0;
}
