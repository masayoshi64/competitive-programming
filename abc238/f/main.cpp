/* #region header */

#pragma GCC optimize("Ofast")
#include <bits/stdc++.h>
using namespace std;
// types
using ll = long long;
using ull = unsigned long long;
using ld = long double;
typedef pair<ll, ll> Pl;
typedef pair<int, int> Pi;
typedef vector<ll> vl;
typedef vector<int> vi;
typedef vector<char> vc;
template <typename T> using mat = vector<vector<T>>;
typedef vector<vector<int>> vvi;
typedef vector<vector<long long>> vvl;
typedef vector<vector<char>> vvc;
// abreviations
#define all(x) (x).begin(), (x).end()
#define rall(x) (x).rbegin(), (x).rend()
#define rep_(i, a_, b_, a, b, ...) for (ll i = (a), max_i = (b); i < max_i; i++)
#define rep(i, ...) rep_(i, __VA_ARGS__, __VA_ARGS__, 0, __VA_ARGS__)
#define rrep_(i, a_, b_, a, b, ...) for (ll i = (b - 1), min_i = (a); i >= min_i; i--)
#define rrep(i, ...) rrep_(i, __VA_ARGS__, __VA_ARGS__, 0, __VA_ARGS__)
#define srep(i, a, b, c) for (ll i = (a), max_i = (b); i < max_i; i += c)
#define SZ(x) ((int)(x).size())
#define pb(x) push_back(x)
#define eb(x) emplace_back(x)
#define mp make_pair
//入出力
#define print(x) cout << x << endl
template <class T> ostream &operator<<(ostream &os, const vector<T> &v) {
    for (auto &e : v)
        cout << e << " ";
    cout << endl;
    return os;
}
void scan(int &a) {
    cin >> a;
}
void scan(long long &a) {
    cin >> a;
}
void scan(char &a) {
    cin >> a;
}
void scan(double &a) {
    cin >> a;
}
void scan(string &a) {
    cin >> a;
}
template <class T> void scan(vector<T> &a) {
    for (auto &i : a)
        scan(i);
}
#define vsum(x) accumulate(all(x), 0LL)
#define vmax(a) *max_element(all(a))
#define vmin(a) *min_element(all(a))
#define lb(c, x) distance((c).begin(), lower_bound(all(c), (x)))
#define ub(c, x) distance((c).begin(), upper_bound(all(c), (x)))
// functions
// gcd(0, x) fails.
ll gcd(ll a, ll b) {
    return b ? gcd(b, a % b) : a;
}
ll lcm(ll a, ll b) {
    return a / gcd(a, b) * b;
}
ll safemod(ll a, ll b) {
    return (a % b + b) % b;
}
template <class T> bool chmax(T &a, const T &b) {
    if (a < b) {
        a = b;
        return 1;
    }
    return 0;
}
template <class T> bool chmin(T &a, const T &b) {
    if (b < a) {
        a = b;
        return 1;
    }
    return 0;
}
template <typename T> T mypow(T x, ll n) {
    T ret = 1;
    while (n > 0) {
        if (n & 1)
            (ret *= x);
        (x *= x);
        n >>= 1;
    }
    return ret;
}
ll modpow(ll x, ll n, const ll mod) {
    ll ret = 1;
    while (n > 0) {
        if (n & 1)
            (ret *= x);
        (x *= x);
        n >>= 1;
        x %= mod;
        ret %= mod;
    }
    return ret;
}

uint64_t my_rand(void) {
    static uint64_t x = 88172645463325252ULL;
    x = x ^ (x << 13);
    x = x ^ (x >> 7);
    return x = x ^ (x << 17);
}
int popcnt(ull x) {
    return __builtin_popcountll(x);
}
template <typename T> vector<int> IOTA(vector<T> a) {
    int n = a.size();
    vector<int> id(n);
    iota(all(id), 0);
    sort(all(id), [&](int i, int j) { return a[i] < a[j]; });
    return id;
}
struct Timer {
    clock_t start_time;
    void start() {
        start_time = clock();
    }
    int lap() {
        // return x ms.
        return (clock() - start_time) * 1000 / CLOCKS_PER_SEC;
    }
};
template <int Mod> struct modint {
    int x;

    modint() : x(0) {
    }

    modint(long long y) : x(y >= 0 ? y % Mod : (Mod - (-y) % Mod) % Mod) {
    }

    modint &operator+=(const modint &p) {
        if ((x += p.x) >= Mod)
            x -= Mod;
        return *this;
    }

    modint &operator-=(const modint &p) {
        if ((x += Mod - p.x) >= Mod)
            x -= Mod;
        return *this;
    }

    modint &operator*=(const modint &p) {
        x = (int)(1LL * x * p.x % Mod);
        return *this;
    }

    modint &operator/=(const modint &p) {
        *this *= p.inverse();
        return *this;
    }

    modint operator-() const {
        return modint(-x);
    }

    modint operator+(const modint &p) const {
        return modint(*this) += p;
    }

    modint operator-(const modint &p) const {
        return modint(*this) -= p;
    }

    modint operator*(const modint &p) const {
        return modint(*this) *= p;
    }

    modint operator/(const modint &p) const {
        return modint(*this) /= p;
    }

    bool operator==(const modint &p) const {
        return x == p.x;
    }

    bool operator!=(const modint &p) const {
        return x != p.x;
    }

    modint inverse() const {
        int a = x, b = Mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            swap(a -= t * b, b);
            swap(u -= t * v, v);
        }
        return modint(u);
    }

    modint pow(int64_t n) const {
        modint ret(1), mul(x);
        while (n > 0) {
            if (n & 1)
                ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    friend ostream &operator<<(ostream &os, const modint &p) {
        return os << p.x;
    }

    friend istream &operator>>(istream &is, modint &a) {
        long long t;
        is >> t;
        a = modint<Mod>(t);
        return (is);
    }

    static int get_mod() {
        return Mod;
    }

    constexpr int get() const {
        return x;
    }
};
template <typename T = int> struct Edge {
    int from, to;
    T cost;
    int idx;

    Edge() = default;

    Edge(int from, int to, T cost = 1, int idx = -1) : from(from), to(to), cost(cost), idx(idx) {
    }

    operator int() const {
        return to;
    }
};

template <typename T = int> struct Graph {
    vector<vector<Edge<T>>> g;
    int es;

    Graph() = default;

    explicit Graph(int n) : g(n), es(0) {
    }

    size_t size() const {
        return g.size();
    }

    void add_directed_edge(int from, int to, T cost = 1) {
        g[from].emplace_back(from, to, cost, es++);
    }

    void add_edge(int from, int to, T cost = 1) {
        g[from].emplace_back(from, to, cost, es);
        g[to].emplace_back(to, from, cost, es++);
    }

    void read(int M, int padding = -1, bool weighted = false, bool directed = false) {
        for (int i = 0; i < M; i++) {
            int a, b;
            cin >> a >> b;
            a += padding;
            b += padding;
            T c = T(1);
            if (weighted)
                cin >> c;
            if (directed)
                add_directed_edge(a, b, c);
            else
                add_edge(a, b, c);
        }
    }
};

/* #endregion*/
// constant
#define inf 1000000000ll
#define INF 4000000004000000000LL
#define mod 998244353ll
#define endl '\n'
const long double eps = 0.000000000000001;
const long double PI = 3.141592653589793;
using mint = modint<mod>;
using vmint = vector<mint>;
/**
 * @brief Formal Power Series
 * @see https://ei1333.github.io/library/math/fps/formal-power-series.cpp
 * @arg modint<mod>, fft (has a method 'multiply')
 * @docs docs/FormalPowerSeries
 */
struct MULT {};
template <typename T, auto &fft> struct FormalPowerSeries : vector<T> {
    using vector<T>::vector;
    using P = FormalPowerSeries;

    P pre(int deg) const {
        return P(begin(*this), begin(*this) + min((int)this->size(), deg));
    }

    P rev(int deg = -1) const {
        P ret(*this);
        if (deg != -1)
            ret.resize(deg, T(0));
        reverse(begin(ret), end(ret));
        return ret;
    }

    void shrink() {
        while (this->size() && this->back() == T(0))
            this->pop_back();
    }

    P operator+(const P &r) const {
        return P(*this) += r;
    }

    P operator+(const T &v) const {
        return P(*this) += v;
    }

    P operator-(const P &r) const {
        return P(*this) -= r;
    }

    P operator-(const T &v) const {
        return P(*this) -= v;
    }

    P operator*(const P &r) const {
        return P(*this) *= r;
    }

    P operator*(const T &v) const {
        return P(*this) *= v;
    }

    P operator/(const P &r) const {
        return P(*this) /= r;
    }

    P operator%(const P &r) const {
        return P(*this) %= r;
    }

    P &operator+=(const P &r) {
        if (r.size() > this->size())
            this->resize(r.size());
        for (int i = 0; i < r.size(); i++)
            (*this)[i] += r[i];
        return *this;
    }

    P &operator-=(const P &r) {
        if (r.size() > this->size())
            this->resize(r.size());
        for (int i = 0; i < r.size(); i++)
            (*this)[i] -= r[i];
        return *this;
    }

    P &operator*=(const P &r) {
        if (this->empty() || r.empty()) {
            this->clear();
            return *this;
        }
        auto ret = fft.multiply(*this, r);
        return *this = {begin(ret), end(ret)};
    }

    P &operator/=(const P &r) {
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int n = this->size() - r.size() + 1;
        return *this = (rev().pre(n) * r.rev().inv(n)).pre(n).rev(n);
    }

    P &operator%=(const P &r) {
        return *this -= *this / r * r;
    }

    pair<P, P> div_mod(const P &r) {
        P q = *this / r;
        return make_pair(q, *this - q * r);
    }

    P operator-() const {
        P ret(this->size());
        for (int i = 0; i < this->size(); i++)
            ret[i] = -(*this)[i];
        return ret;
    }

    P &operator+=(const T &r) {
        if (this->empty())
            this->resize(1);
        (*this)[0] += r;
        return *this;
    }

    P &operator-=(const T &r) {
        if (this->empty())
            this->resize(1);
        (*this)[0] -= r;
        return *this;
    }

    P &operator*=(const T &v) {
        for (int i = 0; i < this->size(); i++)
            (*this)[i] *= v;
        return *this;
    }

    P dot(P r) const {
        P ret(min(this->size(), r.size()));
        for (int i = 0; i < ret.size(); i++)
            ret[i] = (*this)[i] * r[i];
        return ret;
    }

    P operator>>(int sz) const {
        if (this->size() <= sz)
            return {};
        P ret(*this);
        ret.erase(ret.begin(), ret.begin() + sz);
        return ret;
    }

    P operator<<(int sz) const {
        P ret(*this);
        ret.insert(ret.begin(), sz, T(0));
        return ret;
    }

    T operator()(T x) const {
        T r = 0, w = 1;
        for (auto &v : *this) {
            r += w * v;
            w *= x;
        }
        return r;
    }

    P diff() const {
        const int n = (int)this->size();
        P ret(max(0, n - 1));
        for (int i = 1; i < n; i++)
            ret[i - 1] = (*this)[i] * T(i);
        return ret;
    }

    P integral() const {
        const int n = (int)this->size();
        P ret(n + 1);
        ret[0] = T(0);
        for (int i = 0; i < n; i++)
            ret[i + 1] = (*this)[i] / T(i + 1);
        return ret;
    }

    // F(0) must not be 0
    P inv(int deg = -1) const {
        assert(((*this)[0]) != T(0));
        const int n = (int)this->size();
        if (deg == -1)
            deg = n;
        P ret({T(1) / (*this)[0]});
        for (int i = 1; i < deg; i <<= 1) {
            ret = (ret + ret - ret * ret * pre(i << 1)).pre(i << 1);
        }
        return ret.pre(deg);
    }

    // F(0) must be 1
    P log(int deg = -1) const {
        assert((*this)[0] == T(1));
        const int n = (int)this->size();
        if (deg == -1)
            deg = n;
        return (this->diff() * this->inv(deg)).pre(deg - 1).integral();
    }

    P sqrt(
        int deg = -1, const function<T(T)> &get_sqrt = [](T) { return T(1); }) const {
        const int n = (int)this->size();
        if (deg == -1)
            deg = n;
        if ((*this)[0] == T(0)) {
            for (int i = 1; i < n; i++) {
                if ((*this)[i] != T(0)) {
                    if (i & 1)
                        return {};
                    if (deg - i / 2 <= 0)
                        break;
                    auto ret = (*this >> i).sqrt(deg - i / 2, get_sqrt);
                    if (ret.empty())
                        return {};
                    ret = ret << (i / 2);
                    if (ret.size() < deg)
                        ret.resize(deg, T(0));
                    return ret;
                }
            }
            return P(deg, 0);
        }
        auto sqr = T(get_sqrt((*this)[0]));
        if (sqr * sqr != (*this)[0])
            return {};
        P ret{sqr};
        T inv2 = T(1) / T(2);
        for (int i = 1; i < deg; i <<= 1) {
            ret = (ret + pre(i << 1) * ret.inv(i << 1)) * inv2;
        }
        return ret.pre(deg);
    }

    P sqrt(const function<T(T)> &get_sqrt, int deg = -1) const {
        return sqrt(deg, get_sqrt);
    }

    // F(0) must be 0
    P exp(int deg = -1) const {
        if (deg == -1)
            deg = this->size();
        assert((*this)[0] == T(0));
        const int n = (int)this->size();
        if (deg == -1)
            deg = n;
        P ret({T(1)});
        for (int i = 1; i < deg; i <<= 1) {
            ret = (ret * (pre(i << 1) + T(1) - ret.log(i << 1))).pre(i << 1);
        }
        return ret.pre(deg);
    }

    P pow(int64_t k, int deg = -1) const {
        const int n = (int)this->size();
        if (deg == -1)
            deg = n;
        for (int i = 0; i < n; i++) {
            if ((*this)[i] != T(0)) {
                T rev = T(1) / (*this)[i];
                P ret = (((*this * rev) >> i).log() * k).exp() * ((*this)[i].pow(k));
                if (i * k > deg)
                    return P(deg, T(0));
                ret = (ret << (i * k)).pre(deg);
                if (ret.size() < deg)
                    ret.resize(deg, T(0));
                return ret;
            }
        }
        return *this;
    }

    P mod_pow(int64_t k, P g) const {
        P modinv = g.rev().inv();
        auto get_div = [&](P base) {
            if (base.size() < g.size()) {
                base.clear();
                return base;
            }
            int n = base.size() - g.size() + 1;
            return (base.rev().pre(n) * modinv.pre(n)).pre(n).rev(n);
        };
        P x(*this), ret{1};
        while (k > 0) {
            if (k & 1) {
                ret *= x;
                ret -= get_div(ret) * g;
                ret.shrink();
            }
            x *= x;
            x -= get_div(x) * g;
            x.shrink();
            k >>= 1;
        }
        return ret;
    }

    P taylor_shift(T c) const {
        int n = (int)this->size();
        vector<T> fact(n), rfact(n);
        fact[0] = rfact[0] = T(1);
        for (int i = 1; i < n; i++)
            fact[i] = fact[i - 1] * T(i);
        rfact[n - 1] = T(1) / fact[n - 1];
        for (int i = n - 1; i > 1; i--)
            rfact[i - 1] = rfact[i] * T(i);
        P p(*this);
        for (int i = 0; i < n; i++)
            p[i] *= fact[i];
        p = p.rev();
        P bs(n, T(1));
        for (int i = 1; i < n; i++)
            bs[i] = bs[i - 1] * c * rfact[i] * fact[i - 1];
        p = (p * bs).pre(n);
        p = p.rev();
        for (int i = 0; i < n; i++)
            p[i] *= rfact[i];
        return p;
    }
};
#pragma once
/**
 * @brief Number Theoretic Transform
 * @docs docs/NTT.md
 * @param modint
 */
template <typename Mint> struct NTT {
  private:
    vector<Mint> root_pow, root_pow_inv;
    int max_base;
    Mint root; //原始根

    void ntt(vector<Mint> &a) {
        const int n = a.size();
        assert((n & (n - 1)) == 0);
        assert(__builtin_ctz(n) <= max_base);
        for (int m = n / 2; m >= 1; m >>= 1) {
            Mint w = 1;
            for (int s = 0, k = 0; s < n; s += 2 * m) {
                for (int i = s, j = s + m; i < s + m; ++i, ++j) {
                    auto x = a[i], y = a[j] * w;
                    a[i] = x + y, a[j] = x - y;
                }
                w *= root_pow[__builtin_ctz(++k)];
            }
        }
    }

    void intt(vector<Mint> &a) {
        const int n = a.size();
        assert((n & (n - 1)) == 0);
        assert(__builtin_ctz(n) <= max_base);
        for (int m = 1; m < n; m *= 2) {
            Mint w = 1;
            for (int s = 0, k = 0; s < n; s += 2 * m) {
                for (int i = s, j = s + m; i < s + m; ++i, ++j) {
                    auto x = a[i], y = a[j];
                    a[i] = x + y, a[j] = (x - y) * w;
                }
                w *= root_pow_inv[__builtin_ctz(++k)];
            }
        }
    }

  public:
    NTT() {
        const unsigned Mod = Mint::get_mod();
        auto tmp = Mod - 1;
        max_base = 0;
        while (tmp % 2 == 0)
            tmp >>= 1, max_base++;
        root = 2;
        while (root.pow((Mod - 1) >> 1) == 1)
            root += 1;
        root_pow.resize(max_base);
        root_pow_inv.resize(max_base);
        for (int i = 0; i < max_base; i++) {
            root_pow[i] = -root.pow((Mod - 1) >> (i + 2));
            root_pow_inv[i] = Mint(1) / root_pow[i];
        }
    }

    /**
     * @brief 畳み込み
     * @param vector<modint<mod>>
     */
    vector<Mint> multiply(vector<Mint> a, vector<Mint> b) {
        const int need = a.size() + b.size() - 1;
        int nbase = 1;
        while ((1 << nbase) < need)
            nbase++;
        int sz = 1 << nbase;
        a.resize(sz, 0);
        b.resize(sz, 0);
        ntt(a);
        ntt(b);
        Mint inv_sz = Mint(1) / sz;
        for (int i = 0; i < sz; i++)
            a[i] *= b[i] * inv_sz;
        intt(a);
        a.resize(need);
        return a;
    }
};
NTT<mint> ntt;
using FPS = FormalPowerSeries<mint, ntt>;
int main() {
    cin.tie(0);
    ios::sync_with_stdio(0);
    cout << setprecision(30) << fixed;
    ll n, k;
    cin >> n >> k;
    vl p(n), q(n);
    scan(p);
    scan(q);
    Graph<ll> g(n);
    vi root(n, 1);
    rep(i, n) {
        rep(j, n) {
            if (i == j)
                continue;
            if (p[i] < p[j] && q[i] < q[j]) {
                g.add_directed_edge(i, j);
                root[j] = 0;
            }
        }
    }
    vl cnt(n);
    auto rec = [&](int v, auto &rec) -> void {
        for (auto e : g.g[v]) {
            rec(e.to, rec);
            cnt[v] += cnt[e.to];
        }
        cnt[v]++;
        return;
    };
    rep(i, n) {
        if (root[i])
            rec(i, rec);
    }
    auto dp = [&](int v, auto &dp) -> FPS {
        FPS f(1);
        f[0] = 1;
        for (auto e : g.g[v]) {
            f *= dp(e.to, dp);
        }
        f.resize(cnt[v] + 1);
        f[cnt[v]] += 1;
        return f;
    };
    FPS f(1);
    f[0] = 1;
    rep(i, n) {
        if (root[i])
            f *= dp(i, dp);
    }
    f.resize(k + 1);
    print(f[k]);
}