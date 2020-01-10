#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cstdlib>
#include<vector>
#include<queue>
#include<ctime>

using namespace std;

const int maxn=810;
const int maxm=1000010;
const int cutoff_times=6000010;
const int r=500010;

struct Edge{
	int u,v,w;
	int cover;
	int l_pos;
}E[maxm];

struct Vertex{
	bool vis,confChange;
	int dscore,age;
}V[maxn];

int n,m,Clock;
bool b[maxn][maxn],Anti[maxn];
vector<int> G[maxn],C,Ans;
vector<int> L;

void Save_Ans(){
	Ans.clear();
	for(int i=0;i<C.size();++i)
		Ans.push_back(C[i]);
}

void Get_Dscore(){
	
	for(int i=1;i<=m;i++){
		
		if(E[i].cover==1){
			if(V[E[i].u].vis) V[E[i].u].dscore-=E[i].w;
			if(V[E[i].v].vis) V[E[i].v].dscore-=E[i].w;
		}
		else if(E[i].cover==0){
			V[E[i].u].dscore+=E[i].w;
			V[E[i].v].dscore+=E[i].w;
		}
	}
}

void Delete_Max(bool var){
	
	int Max_Dscore=V[C[0]].dscore,Maxi=0;
	for(int i=1;i<C.size();++i)
		if(V[C[i]].dscore>Max_Dscore)
			Max_Dscore=V[C[i]].dscore,Maxi=i;
	
	int u=C[Maxi];
	C[Maxi]=C[C.size()-1],C.pop_back();
	V[u].vis=false;
	
	for(int i=0;i<G[u].size();++i){
		
		int t=G[u][i]; E[t].cover--;
		
		int v=(E[t].u==u)?E[t].v:E[t].u;
		
		if(E[t].cover==1) V[v].dscore-=E[t].w;
		
		else if(E[t].cover==0){
			V[u].dscore+=E[t].w;
			V[v].dscore+=E[t].w;
			
			L.push_back(t);
			E[t].l_pos=L.size()-1; 
		}
	}
	
	if(var){
		V[u].confChange=false;
		for(int i=0;i<G[u].size();i++){
			int t=G[u][i];
			int v=(E[t].u==u)?E[t].v:E[t].u;
			V[v].confChange=true;
		}
	}
}

void Add_Node(int u){
	
	C.push_back(u);
	V[u].vis=true; V[u].age=++Clock; V[u].dscore=0;
	
	for(int i=0;i<G[u].size();++i){
		int t=G[u][i]; E[t].cover++;
		int v=(E[t].u==u)?E[t].v:E[t].u;
		
		if(E[t].cover==2) V[v].dscore+=E[t].w;
		
		else if(E[t].cover==1){
			V[u].dscore-=E[t].w;
			V[v].dscore-=E[t].w;
			
			int pos=E[t].l_pos,temp=L[L.size()-1];
			E[temp].l_pos=pos;
			L[pos]=temp;
			L.pop_back();
		}
	}
}

void Update_w(){
	for(int i=0;i<L.size();++i){
		int t=L[i],oldw=E[t].w;
		E[t].w++;
		
		if(E[t].w>r) E[t].w>>=1;
		
		V[E[t].u].dscore+=E[t].w-oldw;
		V[E[t].v].dscore+=E[t].w-oldw;
	}
}

void Init(){
	
	m=0; Clock=0;
	memset(V,0,sizeof(V));
	Ans.clear(); C.clear(); L.clear();
	for(int i=1;i<=n;i++) G[i].clear();
	
	for(int u=1;u<=n;u++)
		for(int v=u+1;v<=n;v++)
			if(!b[u][v]){
				E[++m]=(Edge){u,v,1,0,0};
				G[u].push_back(m);
				G[v].push_back(m); 
			}
	
	for(int i=1;i<=n;i++) V[i].confChange=true;
	for(int i=m;i>=1;i--){
		int u=E[i].u,v=E[i].v;
		if(!V[u].vis && !V[v].vis){
			
			if(G[u].size()>G[v].size()){
				V[u].vis=true;
				C.push_back(u);  
			}
			else{
				V[v].vis=true;
				C.push_back(v); 
			}
		}
	}
	for(int i=1;i<=m;i++)
		E[i].cover=V[E[i].u].vis+V[E[i].v].vis;
	Get_Dscore();
}

void NuMVC(int cutoff){
	for(int Clock=0;Clock<cutoff;Clock++){
		if(C.size()==0) break;
		 
		if(L.empty()){
			Save_Ans();
			Delete_Max(0);
			continue;
		}
		Delete_Max(1);
		
		int t=L[rand()%(L.size())];
		
		int u=E[t].u,v=E[t].v,p;
		if(V[u].confChange!=V[v].confChange)
			p=V[u].confChange?u:v;
		else if(V[u].dscore!=V[v].dscore)
			p=(V[u].dscore>V[v].dscore)?u:v;
		else
			p=V[u].age<V[v].age?u:v;
		Add_Node(p);
		Update_w();
	}
}

int Read(){
	int x=0;char ch=getchar();
	while(ch>'9' || ch<'0') ch=getchar();
	while(ch>='0' && ch<='9') x=x*10+ch-'0',ch=getchar();
	return x;
}

int main(){
	int u,v;
	srand(time(0));
	while(~scanf("%d%d",&n,&m)){
		
		memset(b,0,sizeof(b));
		for(int i=1;i<=m;i++){
			u=Read();v=Read();
			b[u][v]=b[v][u]=true;
		}
		
		Init(); 
		NuMVC(cutoff_times);
		
		memset(Anti,0,sizeof(Anti));
		for(int i=0;i<Ans.size();++i)
			Anti[Ans[i]]=true;
		
		printf("%d\n",n-Ans.size());
		for(int i=1;i<=n;i++)
			if(!Anti[i])
				printf("%d ",i);
		printf("\n"); 
	}
	return 0;
}
