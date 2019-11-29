#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstring>
using namespace std;
const int MAXN = 805, MAXL = 26;

unsigned int pos[MAXN][2], p[32], n, nlog, curmax;

inline int ones_number(unsigned int n) {
	int ans = 0;
	while (n) {
		n = n & n - 1;
		++ans;
	}
	return ans;
}

struct Bitset {
	unsigned int num[MAXL];
	Bitset(){memset(num, 0, sizeof num);}
	void init(bool a[]) {
		for (int i = 0; i < n; ++i) {
			//printf(a[i]?"1 ":"0 ");
			if (a[i]) {
				num[pos[i][0]] |= p[pos[i][1]];
				//printf("%d %d %d %u\n", i, pos[i][0], pos[i][1], p[pos[i][1]]);
			}
		}
		//puts("");
		//printf("%u\n", num[0]);
	}
	bool is_zeros() {
		for (int i = 0; i < nlog; ++i) if (num[i]) return false;
		return true;
	}
	int get_one_number() {
		int sum = 0;
		for (int i = 0; i < nlog; ++i) sum += ones_number(num[i]);
		return sum;
	}
	int nextbit() {
		for (int i = 0, k = 0; i < nlog; ++i, k += 32)
			if (num[i])
				for (int t = num[i]; t; t >>= 1, ++k)
					if (t&1 == 1)return k;
		return -1;
	}
	void set(int a, bool flag) {
		// printf("%d %d %u ", flag, a, num[pos[a][0]]);
		if (flag) num[pos[a][0]] |= p[pos[a][1]];
		else num[pos[a][0]] &= ~p[pos[a][1]];
		// printf("%u\n", num[pos[a][0]]);
	}
	void operator &= (const Bitset &a) {
		for (int i = 0; i < nlog; ++i) num[i] &= a.num[i];
	}
	Bitset operator & (const Bitset &a) const {
		Bitset ans = *this;
		ans &= a;
		return ans;
	}
	void operator |= (const Bitset &a) {
		for (int i = 0; i < nlog; ++i) num[i] |= a.num[i];
	}
	Bitset operator | (const Bitset &a) const {
		Bitset ans = *this;
		ans |= a;
		return ans;
	}
	Bitset Inv() {
		Bitset ans = *this;
		for (int i = 0; i < nlog; ++i) ans.num[i] = ~ans.num[i];
		return ans;
	}
}neighbor[MAXN], invn[MAXN], ans;

int q_bb[MAXN], c_k[MAXN];
void BBColor(Bitset P, int u[], int c[]) {
	Bitset Q;
	for (int k = 0, col = 1, v; !P.is_zeros(); ++col) {
		Q = P;
		while (!Q.is_zeros()) {
			v = Q.nextbit();
			//printf("%d***\n", v);
			P.set(v, false);
			Q.set(v, false);
			Q &= invn[v];
			u[k] = v;
			c[k++] = col;
		}
	}
}

void save(Bitset C) {
	ans = C;
	curmax = ans.get_one_number();
}

void BBClique(Bitset C, Bitset P) {
	//printf("%d\n", C.get_one_number());
	int m = P.get_one_number();
	int u[m], c[m];
	BBColor(P, u, c);
	for (int i = m - 1, v; ~i; --i) {
		if (c[i] + C.get_one_number() <= curmax) return;
		Bitset Q = P;
		v = u[i];
		//printf("%d ******\n", v);
		C.set(v, true);
		Q &= neighbor[v];
		if (Q.is_zeros() && C.get_one_number() > curmax) save(C);
		if (!Q.is_zeros()) BBClique(C, Q);
		P.set(v, false); C.set(v, false);
	}
}

void init() {
	int cnt0 = 0, cnt1 = 0;
	for (int i = 0; i < n; ++i) {
		pos[i][0] = cnt0, pos[i][1] = cnt1;
		++cnt1;
		if (cnt1 == 32) {
			++cnt0;
			cnt1 = 0;
		}
	}
	p[0] = 1;
	for (int i = 1; i < 32; ++i) p[i] = p[i - 1] << 1;
}

bool mat[MAXN][MAXN], mat2[MAXN][MAXN];
int d[MAXN], index[MAXN], index2[MAXN];

bool cmp(int v, int u) { return d[v] < d[u];}

int main() {
	puts("???");
	freopen("test.clq", "r", stdin);
	int m;
	scanf("%d%d", &n, &m);
	init();
	nlog = (n - 1) / 32 + 1;
	
	for (int i = 1; i < n; ++i) index[i] = i;
	while (m--) {
		static int u, v;
		scanf("%d%d", &u, &v);
		if (u > n || v > n) continue;
		--u, --v;
		++d[u], ++d[v];
		mat2[u][v] = mat2[v][u] = true;
	}
	
	sort(index, index + n, cmp);
	for (int i = 0 ; i < n; ++i) {
		index2[index[i]] = i;
	}
	
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j)
			{
				static int u, v;
				u = index2[i], v = index2[j];
				mat[u][v] = mat[v][u] = mat2[i][j];
			}
	}
	
	for (int i = 0; i < n; ++i) {
		neighbor[i].init(mat[i]);
		for (int j = 0; j < n; ++j) mat[i][j] = !mat[i][j];
		mat[i][i] = false;
		invn[i].init(mat[i]);
	}
	//for (int i = 0; i < n; ++i)
	//	printf("%u\n", neighbor[i].num[0]);
	//return 0;
	
	Bitset C, P;
	bool tmp[n];
	for (int i = 0; i < n; ++i) tmp[i] = true;
	P.init(tmp);
	
	puts("begin");
	BBClique(C, P);
	printf("%d\n", curmax);
	return 0;
}