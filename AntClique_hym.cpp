#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
typedef struct{
    int mode;       // 'c' if strategy is clique; 'v' otherwise
    float* tauV;    // pheromone matrix used when mode='v'
    float** tauE;   // pheromone matrix used when mode='c'
    int nbVertices; // number of vertices in the graph
    int nbEdges;    // number of edges in the graph
    int* degree;    // forall i in 0..nbVertices-1, degree[i] = degree of ith vertex
    int** succ;     // forall i in 0..nbVertices-1, forall j in 0..degree[i]-1, succ[i][j] = jth successor of ith vertex
    int** nonSucc;  // forall i in 0..nbVertices-1, forall j in 0..nbVertices-degree[i]-1, nonsucc[i][j] = jth non successor of ith vertex
} graph;
int iseed;

double randFloat(void) {
    const int ia=16807,ic=2147483647,iq=127773,ir=2836;
    int il,ih,it;
    double rc;
    ih = iseed/iq;
    il = iseed%iq;
    it = ia*il-ir*ih;
    if (it > 0) iseed = it;
    else iseed = ic+it;
    rc = ic;
    return(iseed/rc);
}

int getNextRand(int limit) {
    return (int)(100000*randFloat()) % limit;
}

void seed(unsigned int seed){
    iseed = seed;
}
void createGraph(graph* G) {
    int i, j, nbSucc, nbNonSucc;
    int cp, src, dest;
    int edgesCounter=0;
    G->degree = (int*)calloc(G->nbVertices,sizeof(int));
    G->tauE =(float**)calloc(G->nbVertices,sizeof(float*)); 
    G->succ = (int**)calloc(G->nbVertices,sizeof(int*));
    G->nonSucc = (int**)calloc(G->nbVertices,sizeof(int*));
    for (i=0; i<G->nbVertices; i++){
        G->degree[i] = 0;
        G->tauE[i] = (float*)calloc(G->nbVertices,sizeof(float));
        for (j=0; j<G->nbVertices; j++) G->tauE[i][j] = 0.0;
    }
    for(i=0;i<(G->nbEdges);++i){
        scanf("%d%d",&src,&dest);
        src--; 
        dest--;
        edgesCounter++;
        G->tauE[src][dest]=1.0;
        G->tauE[dest][src]=1.0;
        G->degree[src]++;
        G->degree[dest]++;
    }

    if (edgesCounter != G->nbEdges){
        printf("-- Abnormal number of edges --\n");
        exit(1);
    }
    for (src=0; src<G->nbVertices; src++){
        G->succ[src] = (int*)calloc(G->degree[src],sizeof(int));
        G->nonSucc[src] = (int*)calloc(G->nbVertices-G->degree[src]-1,sizeof(int));
        nbSucc = 0;
        nbNonSucc = 0;
        dest = 0;
        for (dest=0; dest<G->nbVertices; dest++){
            if (G->tauE[src][dest]>0){ G->succ[src][nbSucc] = dest; nbSucc++; }
            else if (src != dest){ G->nonSucc[src][nbNonSucc] = dest; nbNonSucc++; }
        }
    }
}

int isEdge(int i, int j, graph *G){
    // returns true if (i,j) is an edge of G
    return (G->tauE[i][j] > 0);
};

void initPhero(int m, float tauMax, graph* G){
    // Initialize pheromone trails to tauMax and initialize mode to m
    int i, j;
    G->mode = m;
    if (m=='v'){
        G->tauV = (float*)calloc(G->nbVertices,sizeof(float));
        for (i=0; i<G->nbVertices; i++) G->tauV[i] = tauMax;
    }
    else{
        for (i=0; i<G->nbVertices; i++)
            for (j=0; j<G->degree[i]; j++) G->tauE[i][G->succ[i][j]] = tauMax;
    }
}

void evaporate(float rho, float tauMin, graph *G){
    // Evaporate pheromone trails wrt to persistence rate rho
    // ensure that no pheromone trail is lower than tauMin
    int i, j;
    if (G->mode=='v')
        for (i=0; i<G->nbVertices; i++){
            G->tauV[i] *= rho;
            if (G->tauV[i]<tauMin) 
                G->tauV[i]=tauMin;
        }
    else
        for (i=0; i<G->nbVertices; i++)
            for (j=0; j<G->degree[i]; j++){
                G->tauE[i][G->succ[i][j]] *= rho;
                if (G->tauE[i][G->succ[i][j]]<tauMin) 
                    G->tauE[i][G->succ[i][j]] = tauMin;
            }
}


void reinforce(int* clique, int cliqueSize, float qty, float tauMax, graph *G){
    // increase pheromone components associated with clique[0..cliqueSize-1] by qty
    // ensure that no pheromone trail is greater than tauMax
    int i,j;
    if (G->mode=='v')
        for (i=0; i<cliqueSize; i++){
            G->tauV[clique[i]] += qty;
            if (G->tauV[clique[i]]>tauMax) G->tauV[clique[i]]=tauMax;
        }
    else
        for (i=0; i<cliqueSize; i++)
            for (j=i+1; j<cliqueSize; j++){
                G->tauE[clique[i]][clique[j]] += qty;
                if (G->tauE[clique[i]][clique[j]]>tauMax) G->tauE[clique[i]][clique[j]]=tauMax;
                G->tauE[clique[j]][clique[i]] = G->tauE[clique[i]][clique[j]];
            }
};


int compFunc(const void *x, const void *y) {
    // comparison function used to sort vertices in a clique
    int pp,qq;
    int t;
    pp = (int)(*(int *)x);
    qq = (int)(*(int *)y);
    if (pp < qq) t = -1;
    else if (pp == qq) t = 0;
    else  t = 1;
    return t;
}

int displayClique(int cliqueSize, int* clique){
    // displays vertices of clique[0..cliqueSize-1] by increasing order
    int i;
    qsort(clique, cliqueSize, sizeof(int), compFunc);
    printf("%d\n",cliqueSize);
    //printf("\t- clique size: %d vertices\n\t- clique list: ",cliqueSize);
    for(i=0 ; i<cliqueSize ; i++) printf("%d ", clique[i]+1);
    printf("\n");
    //printf("\n--------------------\n");
    return(1);
}

int checkClique (int cliqueSize, int* clique, graph* G) {
    // check that clique[0..cliqueSize-1] actually is a clique
    int i, j;
    for (i=0 ; i<cliqueSize ; i++) {
        for (j=i+1 ; j<cliqueSize ; j++) {
            if (!isEdge(clique[i],clique[j],G)) return(0);
        }
    }
    return(1);
}
float myPow(float x, int y){
    // precondition: y >= 0
    // returns x^y in O(log_2(y))
    float p;
    if (y==0) return 1;
    else if (y==1) return x;
    else if (y==2) return x*x;
    else {
        if (y%2==0){p = myPow(x,y/2); return p*p;}
        else{p = myPow(x,(y-1)/2); return p*p*x;};
    };
}

int choose(float *p, int nbCand){
    /* precondition: nbCand > 0 and for i in 0..nbCand-1, p[i] = sum_{j<i} tau[j]^alpha */
    /* returns k with probability (p[k]-p[k-1])/p[nbCand-1] */
    float f = randFloat();
    int left=0;
    int right=nbCand-1;
    int k;
    float total=p[nbCand-1];
    while (left<right){
        k=(left+right+1)/2;
        if (f<p[k-1]/total) right=k-1;
        else if (f>p[k]/total) left=k+1;
        else return k;
    }
    if (left>=0 && left<nbCand) return left;
    else printf("pb choose\n");
    return k;
}

int walkV(int alpha, int* clique, int firstVertex, graph* G) {
    // input: a graph G, an initial vertex firstVertex, and a parameter alpha
    // output: clique built in a greedy randomized way wrt strategy vertex, starting from initial vertex
    
    int i, v;
    int candidates[G->nbVertices];
    int nbCandidates;
    float total;
    float p[G->nbVertices];
    int cliqueSize = 1;
    
    // initializing the clique with firstVertex
    clique[0]=firstVertex;
    // initializing candidates lists and computing p for all candidates
    nbCandidates = G->degree[firstVertex];
    total = 0;
    for (i=0; i<nbCandidates; i++) {
        candidates[i] = G->succ[firstVertex][i];
        p[i] = myPow(G->tauV[candidates[i]], alpha) + total;
        total = p[i];
    }
    while (nbCandidates>0){
        // choice of the next vertex v within candidates, w.r.t. proba p
        i = choose(p,nbCandidates);
        v = candidates[i];
        // adding v to the clique
        clique[cliqueSize] = v;
        cliqueSize++;
        // removing v from the list of candidates
        nbCandidates--;
        candidates[i] = candidates[nbCandidates];
        // filtering the list of candidates and computing p
        i=0;
        total=0;
        while (i<nbCandidates){
            if ( isEdge(v, candidates[i], G) ) {
                // Candidates[i] is still a candidate -> compute p
                p[i] = myPow(G->tauV[candidates[i]], alpha)+total;
                total = p[i];
                i++;
            }
            else {
                // Candidates[i] is no longer a candidate
                nbCandidates--;
                candidates[i]=candidates[nbCandidates];
            }
        }
    }
    return cliqueSize;
}

int walkE(int alpha, int* clique, int firstVertex, graph* G) {
    // input: a graph G, an initial vertex firstVertex, and a parameter alpha
    // output: clique built in a greedy randomized way wrt strategy clique, starting from initial vertex
    // returns the size of clique
    
    int i, v;
    int candidates[G->degree[firstVertex]];
    int nbCandidates;
    float total;
    float p[G->degree[firstVertex]];
    float tauClique[G->nbVertices];
    int cliqueSize = 1;
    
    // initializing the clique with firstVertex
    clique[0]=firstVertex;
    // initializing candidates lists and computing p for all candidates
    nbCandidates = G->degree[firstVertex];
    total = 0;
    for (i=0; i<nbCandidates; i++) {
        candidates[i] = G->succ[firstVertex][i];
        tauClique[candidates[i]] = (G->tauE)[firstVertex][candidates[i]];
        p[i] = myPow(tauClique[candidates[i]], alpha) + total;
        total = p[i];
    }
    while (nbCandidates>0){
        // choice of the next vertex v within candidates, w.r.t. proba p
        i = choose(p,nbCandidates);
        v = candidates[i];
        // adding v to the clique
        clique[cliqueSize] = v;
        cliqueSize++;
        // removing v from the list of candidates
        nbCandidates--;
        candidates[i] = candidates[nbCandidates];
        // filtering the list of candidates and computing p
        i=0;
        total=0;
        while (i<nbCandidates){
            if ( isEdge(v, candidates[i], G) ) {
                // Candidates[i] is still a candidate -> compute p
                tauClique[candidates[i]] += (G->tauE)[v][candidates[i]];
                p[i] = myPow(tauClique[candidates[i]], alpha)+total;
                total = p[i];
                i++;
            }
            else {
                // Candidates[i] is no longer a candidate
                nbCandidates--;
                candidates[i]=candidates[nbCandidates];
            }
        }
    }
    return cliqueSize;
}


int walk(int alpha, int* clique, int firstVertex, graph* G) {
    if (G->mode=='v') return walkV(alpha,clique,firstVertex,G);
    else return walkE(alpha,clique,firstVertex,G);
}


int repair(int* clique, int* cliqueSize, graph* G) {
    // improves clique by greedy local search
    
    int vi, vj, vk, i, j, k, l, stop, found;
    int delta = 0;
    
    i = getNextRand(*cliqueSize);
    stop = i;
    while (1){
        vi = clique[i];
        // looking for (vj,vk) in nonSucc(vi) so that (vj,vk) in E 
        // and vj and vk are connected to all the vertices of the clique
        found = 0;
        for (j=0; ((j<G->nbVertices-G->degree[vi]-1) && (found==0)); j++){
            vj = G->nonSucc[vi][j];
            for (k=j+1; ((k<G->nbVertices-G->degree[vi]-1) && (found==0)); k++){
                vk = G->nonSucc[vi][k];
                if (isEdge(vj,vk,G)){
                    for (l=0; ((l<*cliqueSize) && 
                               ((clique[l]==vi) || 
                                ((isEdge(vj,clique[l],G)) && (isEdge(vk,clique[l],G))))); 
                         l++);
                    if (l==*cliqueSize) found = 1;
                    else if (!isEdge(vj,clique[l],G)) k = G->nbVertices;
                }
            }
        }
        if (found==1){
            // remove vi from the clique, and add vj and vk
            clique[i] = vj;
            clique[*cliqueSize] = vk;
            (*cliqueSize)++;
            delta++;
            stop = i;
        }
        else{
            i++;
            if (i==*cliqueSize) i=0;
            if (i==stop) return delta;
        }
    }
}

int main(){
    int alpha = 1;          // pheromone factor weight
    int best = 10000;       // bound on the size of the maximum clique
    float rho = 0.01;       // persistence rate
    int maxCycles = 5000;   // maximum number of cycles
    int nbAnts = 30;        // number of ants
    float tauMin = 0.01;    // minimum bound on pheromone trails
    float tauMax = 6;       // maximum bound on pheromone trails
    int mustRepair = 1;     // when mustRepair=1, local search is performed on the best clique of each cycle
    int iseed = 3;          // seed of the random number generator
    char mode = 'v';        // 'v' for vertex strategy; 'c' for clique strategy
    clock_t tClock=clock(); // starting time
    seed(iseed);
    int a,b;
    while(~scanf("%d%d",&a,&b)){
        graph G;
        G.nbVertices = a;
        G.nbEdges = b;
        createGraph(&G);
        //printf("Graph: %d vertices and %d edges\n", G.nbVertices, G.nbEdges);
        initPhero(mode,tauMax,&G);
        if (best>G.nbVertices) 
            best=G.nbVertices;
        int bestCliqueSize = 0;
        int bestClique[G.nbVertices];
        int currentCliqueSize;
        int currentClique[G.nbVertices];
        int bestCliqueCycleSize;
        int bestCliqueCycle[G.nbVertices];
        int nbRepair = 0;
        float qty;
        int i, j;
        int nbIter = 0;
        int nbCycles;
        // start the solution process
        for (nbCycles=0; ((nbCycles<maxCycles) && (bestCliqueSize<best)); nbCycles++){
            bestCliqueCycleSize = 0;
            // Each ant computes a clique
            for(i=0 ; ((i<nbAnts) && (bestCliqueSize<best)) ; i++) {
                currentCliqueSize = walk(alpha,currentClique,getNextRand(G.nbVertices),&G);
                nbIter += currentCliqueSize;
                if (currentCliqueSize > bestCliqueCycleSize) {
                    for (j=0; j<currentCliqueSize; j++) 
                        bestCliqueCycle[j] = currentClique[j];
                    bestCliqueCycleSize = currentCliqueSize;
                }
            }
            if ((mustRepair==1) && (bestCliqueCycleSize<best)) 
                // repair the best clique of the cycle
                nbRepair += repair(bestCliqueCycle,&bestCliqueCycleSize,&G);
            // update best clique
            if (bestCliqueCycleSize > bestCliqueSize) {
                /* 
                 if(!checkClique(bestCliqueCycleSize,bestCliqueCycle,&G)){
                    printf(" * WARNING * new clique is not a clique...\n");
                    displayClique(bestCliqueCycleSize,bestCliqueCycle);
                }
                 */
                bestCliqueSize=bestCliqueCycleSize;
                for (j=0; j<bestCliqueSize; j++) 
                    bestClique[j]=bestCliqueCycle[j];
                //printf("Best clique size = %d at cycle %d (%f seconds)\n", bestCliqueSize, nbCycles, (float)(clock()-tClock)/CLOCKS_PER_SEC);
                if (bestCliqueSize >= best){
                    displayClique(bestCliqueSize, bestClique);
                    return 0;
                }
            }
            // Pheromone updating step
            evaporate(1-rho,tauMin,&G);
            qty = 1.0/(float)(bestCliqueSize + 1 - bestCliqueCycleSize);
            reinforce(bestCliqueCycle,bestCliqueCycleSize, qty, tauMax, &G);
        }
        //printf("(%f seconds)\n", (float)(clock()-tClock)/CLOCKS_PER_SEC);
        displayClique(bestCliqueSize, bestClique);
    }
    
    return 0;
}
