#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
typedef struct{
	int mode;		// 'c' if strategy is clique; 'v' otherwise
	float* tauV;	// pheromone matrix used when mode='v'
	float** tauE;	// pheromone matrix used when mode='c'
	int nbVertices;	// number of vertices in the graph
	int nbEdges;	// number of edges in the graph
	int* degree;	// forall i in 0..nbVertices-1, degree[i] = degree of ith vertex
	int** succ;		// forall i in 0..nbVertices-1, forall j in 0..degree[i]-1, succ[i][j] = jth successor of ith vertex
	int** nonSucc;	// forall i in 0..nbVertices-1, forall j in 0..nbVertices-degree[i]-1, nonsucc[i][j] = jth non successor of ith vertex
} graph;

typedef struct{
	int nbVertices;		// number of vertices in the graph
	int* freq;			// forall i in 0..nbVertices-1, freq[i] = number of times 
						// vertex i has been selected since the last reset
	int intersections;	// = sum_{i in 0..nbVertices-1} freq[i]*(freq[i]-1)/2
	int totVertices;	// sum of the sizes of all cliques built since the last reset
	int nbWalks;		// number of cliques built since the last reset
	int nbConflicts;	// number of conflicts in the hashing table 
						// -> gives a lower bound on the number of times a same clique has been recomputed
	int prime;          // size of the first dimension of hashing table
	int sizeOfMem;		// size of the second dimension of the hashing tables
	int **hashingTable;	// forall i in 0..prime-1, forall j in 0..sizeOfMem-1
						// if hashingTable[i][j] != -1, then a walk has computed a clique c such that
						// the first hash key of c if i and the second hash key of c is hashingTable[i][j]
	int *pow3;			// forall i in 0..nbVertices-1, pow3[i] = 3^i (no care of capacity overflow)
	int *perm;			// perm[0..nbVertices-1] is a permutation of 0..nbVertices used to compute the second hash key
} statistics;

int iseed;

double randFloat(void) {
    const int ia=16807,ic=2147483647,iq=127773,ir=2836;
    int il,ih,it;
    double rc;
    ih = iseed/iq;
    il = iseed%iq;
    it = ia*il-ir*ih;
    if (it > 0)	iseed = it;
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
void createGraph(char* file_name, graph* G) {
	char buffer[1024];
	FILE* fd;
	char type;
	char str[1024];
	int cp, src, dest;
	int edgesCounter=0;
	int i, j, nbSucc, nbNonSucc;
	
	// Opening file
	if ( (fd=fopen(file_name, "r"))==NULL){	
		printf("ERROR: Cannot open ascii input file %s", file_name); 
		exit(1);	
	}
	
	// Reading file
	while( (fgets(buffer, 1024, fd))!=NULL) {
		sscanf(buffer, "%c", &type);
		switch (type){
			case 'c': break; // Skip comments
			case 'p': // Problem description
				cp=sscanf(buffer, "%c%s%d%d", &type, str, &(G->nbVertices), &(G->nbEdges));
				if ( (strstr(str, "edge") != NULL) && (cp == 4) ) {
					G->degree = (int*)calloc(G->nbVertices,sizeof(int));
					G->tauE =(float**)calloc(G->nbVertices,sizeof(float*)); 
					G->succ = (int**)calloc(G->nbVertices,sizeof(int*));
					G->nonSucc = (int**)calloc(G->nbVertices,sizeof(int*));
					for (i=0; i<G->nbVertices; i++){
						G->degree[i] = 0;
						G->tauE[i] = (float*)calloc(G->nbVertices,sizeof(float));
						for (j=0; j<G->nbVertices; j++) G->tauE[i][j] = 0.0;
					}
				} 
				else {
					printf("-- Error while reading problem description --\n");
					exit(1);
				}
				break;
			case 'n': break; // Vertex description
			case 'e': // Edge description
				cp=sscanf(buffer, "%c%d%d", &type, &src, &dest);
				src--; 
				dest--;
				if ((cp==3) && (src>=0) && (src<G->nbVertices) && (dest>=0) && (dest<G->nbVertices)){
					edgesCounter++;
					G->tauE[src][dest]=1.0;
					G->tauE[dest][src]=1.0;
					G->degree[src]++;
					G->degree[dest]++;
				} 
				else {
					printf("-- Error while reading edge description #%d: (%d,%d)\n",edgesCounter,src,dest);
					exit(1);
				}
				break;
			default: break;
		}
	}
	fclose(fd);
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
	printf("--------------------\nClique list dump:\n");
	printf("\t- clique size: %d vertices\n\t- clique list: ",cliqueSize);
	for(i=0 ; i<cliqueSize ; i++) printf("%d ", clique[i]+1);
	printf("\n--------------------\n");
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

double getAverageDispersion(int n, int k){
	// returns (1/k) * sum_{i=1}^k{i*C(k,i)*C(n-k,k-i)/C(n,k)}
	// where C(n,k) is the nb of different sets of k elements that can
	// be built from a set of n elements, i.e., C(n,k) = n!/((n-k)!k!)
	int i, j;
	double p_i, prod;
	if (2*k-n < 1) i=1; else i=2*k-n;
	p_i = (double)(k*k)/(double)(i);
	for (j=1; j<i; j++) p_i = p_i*(double)(k-j)*(double)(k-j)/(double)(j);
	for (j=k-i; j<=k-1; j++) p_i = p_i/(double)(n-j);
	for (j=0; j<k-i; j++) p_i = p_i*(double)(n-k-j)/(double)(n-j);
	prod = (double)(i)*p_i;
	//  printf("i=%d, p_i=%f et prod=%f\n", i,p_i,prod);
	for (i++; i<=k; i++){
		p_i = p_i*(double)((k-i+1)*(k-i+1))/(double)(i*(n-2*k+i));
		// p_i = C(k,i)*C(n-k,k-i)/C(n,k)
		//     = probability that |S1 intersection S2| = i 
		//       where S1 and S2 are 2 subsets of {1,2,...,n} such that |S1|=|S2|
		//     = (k^2/i) * pi_{j=1}^{i-1}{(k-j)^2/j} 
		//               * pi_{j=0}{k-i-1}{(n_k-j)/(n-j)}
		//               * pi_{j=k-i}^{k-1}{1/(n-j)}
		prod += (double)(i)*p_i;
		// printf("i=%d, p_i=%f et prod=%f\n", i,p_i,prod);
	}
	prod = prod/(double)(k);
	//  printf("prod = %f\n",prod);
	return prod;
}

void reset(statistics *S){
	// reset data structures used to compute the dispersion rate and the average quality
	int i;
	for (i=0; i<S->nbVertices; i++) S->freq[i]=0;
	S->intersections = 0;
	S->totVertices = 0;
	S->nbWalks = 0;
}

float getSimilarity(statistics *S){
	// return the average similarity of all cliques computed since the last reset
	if (S->nbWalks>1) 
		return (float)(2*S->intersections)/(float)((S->nbWalks-1)*S->totVertices);
	else return 0;
}

float getAverageSize(statistics *S){
	// returns the average size of all cliques since the last reset
	if (S->nbWalks>0) 
		return (float)(S->totVertices)/(float)(S->nbWalks);
	else return 0;
}

int getNbConflicts(statistics *S){
	// returns the number of conflicts that have occur in the hashing table
	// this gives an upper bound on the number of times a same clique has been recomputed
	return S->nbConflicts;
}

void createStatistics(int nbV, statistics *S){
	// initialize data structure used to compute statistics
	int i,j,aux;
	
	S->nbConflicts = 0;
	S->nbVertices = nbV;
	S->freq = (int*)calloc(nbV, sizeof(int));
	reset(S);
	
	S->prime = 452279;
	S->sizeOfMem = 10;
	// initialization of the hashing hashing table: 
	// forall i in 0..prime-1, foralll j in 0..sizeOfMem-1, hashingTable[i][j]=-1
	S->hashingTable = (int**)calloc( S->prime, sizeof(int*) );
	for (i=0; i<S->prime; i++){
		S->hashingTable[i] = (int*)calloc( S->sizeOfMem, sizeof(int) );
		for (j=0; j<S->sizeOfMem; j++) 
			S->hashingTable[i][j]=-1;
	}
	// initialization of pow3: forall i in 0..nbVertices-1, pow3[i]=(2^i) % prime
	S->pow3 = (int*)calloc( (S->nbVertices), sizeof(int) );
	S->pow3[0]= 1;
	for (i=1; i<S->nbVertices; i++) 
		S->pow3[i] = S->pow3[i-1]*3;
	// initialization of perm[0..nbVertices-1] to a permutation of [0..nbVertices-1]
	S->perm = (int*)calloc( S->nbVertices, sizeof(int) );
	for (i=0; i<S->nbVertices; i++) S->perm[i]=i;
	for (i=0; i<S->nbVertices; i++){
		j = i+getNextRand(S->nbVertices-i);
		aux = S->perm[i];
		S->perm[i] = S->perm[j];
		S->perm[j] = aux;
	}
}


void update(int* clique, int cliqueSize, statistics *S){
	// update statistics wrt clique[0..cliqueSize-1]
	int i;
	int h1 = 0;
	int h2 = 0;
	
	S->nbWalks++;
	for (i=0; i<cliqueSize; i++){
		S->intersections += S->freq[clique[i]];
		S->freq[clique[i]]++;
		h1 = h1+S->pow3[clique[i]];
		h2 = h2+S->pow3[S->perm[clique[i]]];
	}
	if (h1 < 0) h1 = -h1;
	h1 = h1 % S->prime;
	S->totVertices += cliqueSize;
	for (i=0; ((i<S->sizeOfMem) && (S->hashingTable[h1][i]!=h2) && (S->hashingTable[h1][i]>=0)); i++);
	if (i==S->sizeOfMem)
		printf("Memory exceeded... should realloc !\n");
	else if (S->hashingTable[h1][i]==h2) 
		S->nbConflicts++; // conflict detected
	else // hashingTable[h1][i]<0) -> add (h1,h2) to the hashing hashingTablele
		S->hashingTable[h1][i]=h2;
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



/**********************/
int usage (char *exec) {
	printf("\nUsage :\n\n\t%s ",exec);
	printf("-a (alpha: int) -B (size of the maximum clique: int) -r (rho: float) ");
	printf("-c (max nb cycle: int) -n (nb ants: int) -m (tau min: float) ");
	printf("-M (tau max: float) -i (DIMACS graph filename) -v display-frequency ");
	printf("[-p (mustRepair set to 1)] -s (seed: positive int) -t (mode: v if AntClique(Vertex); c if AntClique(Clique))\n\n");
	return(0);
}

/**********************/
int parse(int* alpha, 
		  int* best, 
		  float* rho, 
		  int* maxCycles, 
		  int* nbAnts, 
		  float* tauMin, 
		  float* tauMax, 
		  int* verbose, 
		  int* displayFreq, 
		  char* fileName, 
		  int* mustRepair, 
		  int* iseed, 
		  char* mode, 
		  char* argv[], int argc){
	// set parameters wrt argv and argc
	char ch;
	extern char* optarg;
	while ( (ch = getopt(argc, argv, "a:B:r:p:c:n:m:M:v:i:?:h:s:t:"))!=-1 ) {
		switch(ch) {
			case 'a': *alpha=atoi(optarg); break;
			case 'B': *best=atoi(optarg); break;
			case 'r': *rho=atof(optarg); break;
			case 'p': *mustRepair=atoi(optarg); break;
			case 's': *iseed=atoi(optarg); break;
			case 'c': *maxCycles=atoi(optarg); break;
			case 'n': *nbAnts=atoi(optarg); break;
			case 'm': *tauMin=atof(optarg); break;
			case 'M': *tauMax=atof(optarg); break;
			case 'v': *verbose=1; *displayFreq=atoi(optarg); break;
			case 'i': strncpy(fileName, optarg, 254); break;
			case 't': *mode=optarg[0]; break;
			case '?':
			case 'h':
			default: usage(argv[0]); return(1);
		}
	}
	return(0);
}



int main (int argc, char *argv[]) {
	
	// Declare parameters and initialize them with default values
	char fileName[1024];	// name of the file which contains the graph
	int verbose = 0;		// when verbose=1, information is displayed during the solution process
	int displayFreq = 200;  // display frequency when verbose = 1
	int alpha = 1;          // pheromone factor weight
	int best = 10000;		// bound on the size of the maximum clique
	float rho = 0.01;		// persistence rate
	int maxCycles = 3000;	// maximum number of cycles
	int nbAnts = 30;		// number of ants
	float tauMin = 0.01;	// minimum bound on pheromone trails
	float tauMax = 6;		// maximum bound on pheromone trails
	int mustRepair = 1;		// when mustRepair=1, local search is performed on the best clique of each cycle
	int iseed = 3;			// seed of the random number generator
	char mode = 'v';		// 'v' for vertex strategy; 'c' for clique strategy
	
	clock_t tClock=clock();	// starting time
	
	// get parameters
	if( parse(&alpha, 
			  &best, 
			  &rho, 
			  &maxCycles, 
			  &nbAnts, 
			  &tauMin, 
			  &tauMax, 
			  &verbose, 
			  &displayFreq, 
			  fileName, 
			  &mustRepair, 
			  &iseed, 
			  &mode, 
			  argv, argc) == 1) 
		return(1);
	if (verbose==1)
		printf("Params: alpha=%i best=%d rho=%f tauMin=%f tauMax=%f nbCycles=%d nbAnts=%d verbose=%d(%d) input=%s mustRepair=%d seed=%d mode=%c\n",alpha, best, rho, tauMin, tauMax, maxCycles, nbAnts, verbose, displayFreq, fileName, mustRepair, iseed, mode);
	
	// initialize data structures
	seed(iseed);
	graph G;
	createGraph(fileName,&G);
	if (verbose==1) 
		printf("Graph: %d vertices and %d edges\n", G.nbVertices, G.nbEdges);
	initPhero(mode,tauMax,&G);
	if (verbose==1) 
		printf("Pheromone trails initialized to %f\n",tauMax);
	if (best>G.nbVertices) 
		best=G.nbVertices;
	statistics S;
	createStatistics(G.nbVertices,&S);
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
	if (verbose==1) 
		printf("Starting solution process\n");
	for (nbCycles=0; ((nbCycles<maxCycles) && (bestCliqueSize<best)); nbCycles++){
		bestCliqueCycleSize = 0;
		// Each ant computes a clique
		for(i=0 ; ((i<nbAnts) && (bestCliqueSize<best)) ; i++) {
			currentCliqueSize = walk(alpha,currentClique,getNextRand(G.nbVertices),&G);
			if (verbose==1) 
				update(currentClique,currentCliqueSize,&S);
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
			printf("Best clique size = %d at cycle %d (%f seconds)\n", bestCliqueSize, nbCycles, (float)(clock()-tClock)/CLOCKS_PER_SEC);
			if (bestCliqueSize >= best){
				displayClique(bestCliqueSize, bestClique);
				return 0;
			}
		}
		// Pheromone updating step
		evaporate(1-rho,tauMin,&G);
		qty = 1.0/(float)(bestCliqueSize + 1 - bestCliqueCycleSize);
		reinforce(bestCliqueCycle,bestCliqueCycleSize, qty, tauMax, &G);
		// Display
		if((verbose==1) && ((nbCycles+1) % displayFreq == 0)) {
			printf("Statistics at cycle %d : ", nbCycles+1);
			printf("CPU time=%f ",(float)(clock()-tClock)/CLOCKS_PER_SEC);
			printf("Number of LS moves=%d ",nbRepair);
			printf("Average clique size=%f ",getAverageSize(&S));
			printf("Average similarity=%f ",getSimilarity(&S));
			printf("Number of resamplings=%d\n",getNbConflicts(&S));
			reset(&S);
		}
	}
	displayClique(bestCliqueSize, bestClique);
	return(0);
}

