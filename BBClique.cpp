#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
using namespace std;
const int MAXN = 1005, MAXL = 32;
// int MV = 100000000;

clock_t startTime, endTime;
unsigned int pos[MAXN][2], p[32], n, nlog, curmax = 1;

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
        for (int i = 0; i < n; ++i)
            if (a[i])
                num[pos[i][0]] |= p[pos[i][1]];
    }
    inline bool is_zeros() {
        for (int i = 0; i < nlog; ++i) if (num[i]) return false;
        return true;
    }
    inline int get_one_number() {
        int sum = 0;
        for (int i = 0; i < nlog; ++i) sum += ones_number(num[i]);
        return sum;
    }
    inline int nextbit() {
        for (int i = 0, k = 0; i < nlog; ++i, k += 32)
            if (num[i])
                for (int t = num[i]; t; t >>= 1, ++k)
                    if (t&1 == 1)return k;
        return -1;
    }
    inline void set(int a, bool flag) {
        if (flag) num[pos[a][0]] |= p[pos[a][1]];
        else num[pos[a][0]] &= ~p[pos[a][1]];
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
}neighbor[MAXN], invn[MAXN], ans;

inline void outans() {
    printf("%d\n", curmax);
    bool flag = true;
    for (int i = 0; i < n; ++i) {
        if(ans.num[pos[i][0]] & p[pos[i][1]]) {
            if (flag){printf("%d", i + 1);flag = false;}
            else printf(" %d", i + 1);
        }
    }
    exit(0);
}

int tot;
int q_bb[MAXN], c_k[MAXN];
void BBColor(Bitset P, int u[], int c[]) {
    Bitset Q;
    for (int k = 0, col = 1, v; !P.is_zeros(); ++col) {
        // if (tot >= MV) outans();
        Q = P;
        while (!Q.is_zeros()) {
            ++tot;
            v = Q.nextbit();
            P.set(v, false);
            Q.set(v, false);
            Q &= invn[v];
            u[k] = v;
            c[k++] = col;
        }
    }
}

inline void save(Bitset& C) {
    ans = C;
    curmax = ans.get_one_number();
    printf("save %d\n", curmax);
    endTime = clock();
    cout<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<'\n';
}

void BBClique(Bitset C, Bitset P) {
    //if (tot >= MV) outans();
    int m = P.get_one_number();
    int u[m], c[m];
    BBColor(P, u, c);
    for (int i = m - 1, v; ~i; --i) {
        ++tot;
        if (c[i] + C.get_one_number() <= curmax) return;
        Bitset Q = P;
        v = u[i];

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

bool mat[MAXN][MAXN];
int d[MAXN], index[MAXN], index2[MAXN];

bool cmp(int v, int u) { return d[v] < d[u];}

int main() {
    startTime = clock();
    char s[4];
    freopen("data/frb30-15-1.clq", "r", stdin);
    int m;
    scanf("%s", s);
    scanf("%s%d%d", s, &n, &m);
    init();
    nlog = (n - 1) / 32 + 1;
    //MV /= nlog;

    //curmax = 44;
    while (m--) {
        static int u, v;
        scanf("%s%d%d", s, &u, &v);
        --u, --v;
        mat[u][v] = mat[v][u] = true;
    }

    for (int i = 0; i < n; ++i) {
        neighbor[i].init(mat[i]);
        for (int j = 0; j < n; ++j) mat[i][j] = !mat[i][j];
        mat[i][i] = false;
        invn[i].init(mat[i]);
    }

    Bitset C, P;
    bool tmp[n];
    for (int i = 0; i < n; ++i) tmp[i] = true;
    P.init(tmp);

    BBClique(C, P);
    outans();
}
