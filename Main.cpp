/*
* @Author: mkg
* @Date:   2016-05-14 11:54:33
* @Last Modified by:   research1
* @Last Modified time: 2016-09-05 23:40:51
*/

#include <map>
#include <unordered_map>
#include <set>
#include <list>
#include <cmath>
#include <ctime>
#include <deque>
#include <queue>
#include <stack>
#include <bitset>
#include <cstdio>
#include <limits>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <numeric>
#include <utility>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


typedef long long Long;
typedef double DD;
typedef vector<int> VI;
typedef vector<VI > VVI;
typedef pair<int, int> PII;
typedef pair<Long, Long> PLL;
typedef vector<PII> VPII;
typedef vector<PLL> VPLL;

const Long INF = 200000000000000000;
const int MOD = 1000000007;
const Long L_MAX = 9223372036854775807;
const int I_MAX = 2147483647;

#define sf scanf
#define pf printf
#define mem(a,b)          memset(a,b,sizeof(a))
#define pb push_back
#define REP(i,a,b)        for(int i=a; i<=b; ++i)
#define REPI(i,a,b,c)     for(int i=a; i<=b; i+=c)
#define REPR(i,a,b)       for(int i=b; i>=a; --i)
#define REPRD(i,a,b,c)    for(int i=b; i>=a; i-=c)
#define REPB(i,a)         for(int i=a; ;i++)
#define REPRB(i,a)        for(int i=a; ; i--)
#define mp(a,b)   make_pair(a,b)
#define fs        first
#define sc        second
#define SZ(s)     ((int)s.size())
#define PI        3.141592653589793
#define VS        vector<string>
#define VI        vector<int>
#define VD        vector<DD>
#define VL        vector<Long>
#define VVL       vector<VL >
#define lim       10000010
#define tlim      (1<<((int)ceil(log2(lim))+1))
#define unq(vec)  stable_sort(vec.begin(),vec.end());\
                  vec.resize(distance(vec.begin(),unique(vec.begin(),vec.end())));
#define BE(a)     a.begin(),a.end()
#define rev(a)    reverse(BE(a))
#define sorta(a)  stable_sort(BE(a))
#define sortc(a, comp)  sort(BE(a),comp)
#define inf (1<<28)
//int X[]={1,1,2,2,-1,-1,-2,-2},Y[]={2,-2,1,-1,2,-2,1,-1};//knight move
//int X[]={0,-1,-1,-1,0,1,1,1},Y[]={-1,-1,0,1,1,1,0,-1};//8 move
//int X[]={-1,0,1,0},Y[]={0,1,0,-1};//4 move
#define piii pair< int, pair<int,int> >

#include <string.h>
#include <ctype.h>
#include "kseq.h"
#include <zlib.h>


/////// These header files are included for Thread///////////////
#include <chrono>                                              //
#include "ThreadPool.h"                                        //
#include <chrono>                                              //
/////////////////////////////////////////////////////////////////


int K = 14; ///minimizer length
int W = 24; ///minimizer windows size
int ADK;
Long kmask;
Long pat2num(string s)
{
    Long ret = 0;
    int now = 0, kk = 0;
    for (int i = s.size() - 1; i >= 0; i--)
    {
        kk = 0;
        if (s[i] == 'G' || s[i] == 'T') kk |= 2;
        if (s[i] == 'C' || s[i] == 'T') kk |= 1;
        ret |= (kk << now);
        now += 2;

    }
    return ret;
}

Long pat2num2(string s)
{
    Long ret = 0;
    int now = 0, kk = 0;
    for (int i = s.size() - 1; i >= 0; i--)
    {
        if (s[i] == '_') continue;
        kk = 0;
        if (s[i] == 'G' || s[i] == 'T') kk |= 2;
        if (s[i] == 'C' || s[i] == 'T') kk |= 1;
        ret |= (kk << now);
        now += 2;
    }
    return ret;
}

string num2pat(Long num, int k)
{
    string s = "";
    int tem = 0;
    for (int i = 0; i < k; i++)
    {
        tem = num & 3;
        if (tem == 0) s += 'A';
        else if (tem == 1) s += 'C';
        else if (tem == 2) s += 'G';
        else s += 'T';
        num >>= 2;
    }
    reverse(s.begin(), s.end());
    return s;
}

int nt2num(char x)
{
    if (x == 'A') return 0;
    else if (x == 'C') return 1;
    else if (x == 'G') return 2;
    return 3;
}

char num2nt(int x)
{
    if (x & 1)
    {
        if (x & 2) return 'T';
        else return 'C';
    }
    else
    {
        if (x & 2) return 'G';
        else return 'A';
    }
}

//unordered_map<char, char> base;
Long rev_comp(Long x, int k)
{
    Long rev = 0;
    for (int i = 0; i < k; i++)
    {
        rev <<= 2;
        rev |= ( (x & 3) ^ 3);
        x >>= 2;
    }
    return rev;
}

///read,reference
string REFF;

struct my_data
{
    Long minim;
    int idx;
    my_data(Long t, int p): minim(t), idx(p) {}
};

unordered_map<Long, vector< piii > > retf, retr;
vector<int> for_index, rev_index;

void find_minimizers_forward(string s)
{
    retf.clear();
    for_index.clear();
    if (SZ(s) < K or SZ(s) < W)
    {
        if (SZ(s) < W) cerr << "Sequence Length is smaller than W after inserting gap(-).\n";
        else cerr << "Sequence Length is smaller than K after inserting gap(-).\n";
        return ;
    }
    int n = SZ(s), pos = -1, st, en;
    Long temps, minimizer, prevmin = -1;
    piii toinsert;
    vector<piii > v;
    deque<my_data> sliding_window;
    temps = pat2num2(s.substr(0, K));
    sliding_window.push_back(my_data(temps, 0));
    int i ;
    //cerr << num2pat(temps, ((K / 3) * 2)) << " " << temps << endl;
    for (i = K; i < n && i < W; ++i)
    {
        if (s[i] != '_')
        {
            temps <<= 2;
            temps &= kmask;
            temps |= nt2num(s[i]);
        }
        minimizer = temps;
        //cerr << num2pat(temps, ((K / 3) * 2)) << "++" << num2pat(minimizer, ((K / 3) * 2)) << endl;
        while ((!sliding_window.empty()) && minimizer < sliding_window.back().minim)
            sliding_window.pop_back();
        sliding_window.push_back(my_data(minimizer, (i - K + 1)));
    }
    st = 0;
    en = i - 1;
    prevmin = sliding_window.front().minim;
    pos = sliding_window.front().idx;
    //cerr << i << "++" << endl;
    for (; i < n; i++)
    {
        if (s[i] != '_')
        {
            temps <<= 2;
            temps &= kmask;
            temps |= nt2num(s[i]);
        }
        if (s[i - K] != '_') continue;
        minimizer = temps;
        while ((!sliding_window.empty()) && sliding_window.front().idx <= (i - W))
        {

            sliding_window.pop_front();
        }
        while ((!sliding_window.empty()) && minimizer < sliding_window.back().minim)
            sliding_window.pop_back();
        sliding_window.push_back(my_data(minimizer, (i - K + 1)));

        if (sliding_window.front().minim == prevmin)
        {
            en = i - 1;
            pos = sliding_window.front().idx;
        }
        else
        {
            toinsert.first = pos;
            toinsert.second.first = st;
            toinsert.second.second = en;
            if (retf.count(prevmin) == 0)
            {
                v.clear();
                v.emplace_back(toinsert);
                retf[ prevmin ] = v;
            }
            else retf[ prevmin ].emplace_back(toinsert);
            prevmin = sliding_window.front().minim;
            st = i - W + 1;
            en = i;
            pos = sliding_window.front().idx;
        }

    }
    if (sliding_window.front().minim == prevmin)
    {
        en = i - 1;
        pos = sliding_window.front().idx;
    }
    toinsert.first = pos;
    toinsert.second.first = st;
    toinsert.second.second = en;
    if (retf.count(sliding_window.front().minim) == 0)
    {
        v.clear();
        v.emplace_back(toinsert);
        retf[ prevmin ] = v;
    }
    else retf[ sliding_window.front().minim ].emplace_back(toinsert);
    for (auto t : retf)
    {
        for (auto tt : t.sc) for_index.emplace_back(tt.fs);
    }
    sort(BE(for_index));
    return ;
}

void find_minimizers_reverse(string s)
{
    retr.clear();
    rev_index.clear();
    if (SZ(s) < K or SZ(s) < W)
    {
        if (SZ(s) < W) cerr << "Sequence Length is smaller than W after inserting gap(-).\n";
        else cerr << "Sequence Length is smaller than K after inserting gap(-).\n";
        return ;
    }
    int n = SZ(s), pos = -1, st, en;
    Long temps, minimizer, prevmin = -1;
    piii toinsert;
    vector<piii > v;
    deque<my_data> sliding_window;
    temps = pat2num2(s.substr(0, K));
    sliding_window.push_back(my_data(temps, 0));
    int i ;
    //cerr << num2pat(temps, ((K / 3) * 2)) << " " << temps << endl;
    for (i = K; i < n && i < W; ++i)
    {
        if (s[i] != '_')
        {
            temps <<= 2;
            temps &= kmask;
            temps |= nt2num(s[i]);
        }
        minimizer = temps;
        //cerr << num2pat(temps, ((K / 3) * 2)) << "++" << num2pat(minimizer, ((K / 3) * 2)) << endl;
        while ((!sliding_window.empty()) && minimizer < sliding_window.back().minim)
            sliding_window.pop_back();
        sliding_window.push_back(my_data(minimizer, (i - K + 1)));
    }
    st = 0;
    en = i - 1;
    prevmin = sliding_window.front().minim;
    pos = sliding_window.front().idx;
    //cerr << i << "++" << endl;
    for (; i < n; i++)
    {
        if (s[i] != '_')
        {
            temps <<= 2;
            temps &= kmask;
            temps |= nt2num(s[i]);
        }
        if (s[i - K] != '_') continue;
        minimizer = temps;
        while ((!sliding_window.empty()) && sliding_window.front().idx <= (i - W))
        {

            sliding_window.pop_front();
        }
        while ((!sliding_window.empty()) && minimizer < sliding_window.back().minim)
            sliding_window.pop_back();
        sliding_window.push_back(my_data(minimizer, (i - K + 1)));

        if (sliding_window.front().minim == prevmin)
        {
            en = i - 1;
            pos = sliding_window.front().idx;
        }
        else
        {
            toinsert.first = pos;
            toinsert.second.first = st;
            toinsert.second.second = en;
            if (retr.count(prevmin) == 0)
            {
                v.clear();
                v.emplace_back(toinsert);
                retr[ prevmin ] = v;
            }
            else retr[ prevmin ].emplace_back(toinsert);
            prevmin = sliding_window.front().minim;
            st = i - W + 1;
            en = i;
            pos = sliding_window.front().idx;
        }

    }
    if (sliding_window.front().minim == prevmin)
    {
        en = i - 1;
        pos = sliding_window.front().idx;
    }
    toinsert.first = pos;
    toinsert.second.first = st;
    toinsert.second.second = en;
    if (retr.count(sliding_window.front().minim) == 0)
    {
        v.clear();
        v.emplace_back(toinsert);
        retr[ prevmin ] = v;
    }
    else retr[ sliding_window.front().minim ].emplace_back(toinsert);

    for (auto t : retr)
    {
        for (auto tt : t.sc) rev_index.emplace_back(tt.fs);
    }
    sort(BE(rev_index));
    return ;
}

void insert_gap()
{
    int L = SZ(REFF) - 1;
    REP(i, 0, L)
    {
        if ((i + 1) % 3 == 0 && i != 0) REFF[i] = '_';
    }
    return ;
}

void reverse_comp()
{
    int L = SZ(REFF) - 1;
    for (int i = 0; i <= L ; ++i)
    {
        if (REFF[i] == '_') continue;
        if (REFF[i] == 'A')  REFF[i] = 'T';
        else if (REFF[i] == 'T')  REFF[i] = 'A';
        else if (REFF[i] == 'G')  REFF[i] = 'C';
        else if (REFF[i] == 'C')  REFF[i] = 'G';

    }
    reverse(REFF.begin(), REFF.end());
}

string insert_gap_read(string str)
{
    int L = SZ(str) - 1;
    REP(i, 0, L)
    {
        if ((i + 1) % 3 == 0 && i != 0) str[i] = '_';
    }
    return str;
}
int kmer_found, nxt_found;
void find_kmers_of_read_in_refer(int flag, int offset, int start_in_ref, string read_name, string red_seq)
{
    string kkmer;
    int L = SZ(red_seq);
    int end_in_ref = start_in_ref + offset  + L, counter;
    start_in_ref += offset;
    int lidx, ridx;
    if (flag)
    {
        lidx = lower_bound(BE(for_index), start_in_ref) - for_index.begin();
        ridx = upper_bound(BE(for_index), end_in_ref) - for_index.begin();
        //cerr<<lidx<<" and "<<ridx<<" ++ "<<start_in_ref<<" -- "<<end_in_ref<<" :: "<<SZ(for_index)<<endl;
        if (lidx >= ridx) counter = 0;
        else counter = (ridx - lidx + 1);
    }
    else
    {
        lidx = lower_bound(BE(rev_index), start_in_ref) - rev_index.begin();
        ridx = upper_bound(BE(rev_index), end_in_ref) - rev_index.begin();
        //cerr<<lidx<<" or "<<ridx<<" ++ "<<start_in_ref<<" -- "<<end_in_ref<<" :: "<<SZ(rev_index)<<endl;
        if (lidx >= ridx) counter = 0;
        else counter = (ridx - lidx + 1);
    }
    //FILE* ffp = fopen((char*)read_name.c_str(), "w+");
    cout << read_name << endl;
    kmer_found = -1, nxt_found = 0;
    int mxm = 1;
    for (int ii = 0 ; ii <= L - K; ++ii)
    {
        kkmer = red_seq.substr(ii, K);
        kkmer = insert_gap_read(kkmer);
        Long kmer = pat2num2(kkmer);
        //printf("%s %lld\n", num2pat(kmer, ADK).c_str(), kmer);
        if (flag)
        {
            if (retf.find(kmer) != retf.end())
            {
                nxt_found++;
                if (kmer_found == -1)
                {
                    kmer_found = ii;
                }
                else
                {
                    mxm = max(mxm, (ii - kmer_found - 1));
                    kmer_found = ii;
                }
//                for (auto t : retf[kmer])
//                {
//                    printf("%d ", t.fs);
//                }
//                printf("\n");
            }
//            else printf("Not found in reference.\n");
        }
        else
        {
            if (retr.find(kmer) != retr.end())
            {
                nxt_found++;
                if (kmer_found == -1)
                {
                    kmer_found = ii;
                }
                else
                {
                    mxm = max(mxm, (ii - kmer_found - 1));
                    kmer_found = ii;
                }
//                for (auto t : retr[kmer])
//                {
//                    printf("%d ", t.fs);
//                }
//                printf("\n");
            }
//            else printf("Not found in reference.\n");
        }
    }
    cout << start_in_ref << " " << end_in_ref << " " << counter << "\n";
    cout << nxt_found << "\n";

}

void init()
{
    for (int i = 0; i < ADK; i++)
    {
        kmask <<= 2;
        kmask |= 3;
    }
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, const char **argv)
{
    //ios::sync_with_stdio(false);
    //freopen("/Users/mkg/Desktop/Read Mapping Testing/gg_reads.txt", "r", stdin);
    //freopen("output.txt","w",stdout);
    ADK = 10;
    init();
    ThreadPool pool(1);
    vector< future<int> > results;
    struct timespec start, finish;
    double elapsed, total;
    int L;
    total = 0.0;
    string r_str, rr_str;
    if (argc <= 2)
    {
        cerr << "Please! Give valid file path![optional K and W]\n";
        return 0;
    }
    if (argc == 5)
    {
        r_str = string(argv[1]);
        rr_str = string(argv[2]);
        K = atoi(argv[3]);
        W = atoi(argv[4]);
    }
    else if (argc == 4)
    {
        r_str = string(argv[1]);
        rr_str = string(argv[2]);
        K = atoi(argv[3]);
    }
    else if (argc == 3)
    {
        r_str = string(argv[1]);
        rr_str = string(argv[2]);
    }
    string ref_name ;
    int ref_len;
    gzFile fp;
    kseq_t *seq;
    FILE* ttp = fopen((char*)r_str.c_str(), "r");
    fp = gzdopen(fileno(ttp), "r");
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0)
    {
        REFF = (string)seq->seq.s;
        ref_name = (string)seq->name.s;
        ref_len = seq->seq.l;
    }
    kseq_destroy(seq);
    gzclose(fp);
    fclose(ttp);
    cerr << "Seq length: " << ref_len << endl;
    insert_gap();
    clock_gettime(CLOCK_MONOTONIC, &start);
    find_minimizers_forward(REFF);
    reverse_comp();
    find_minimizers_reverse(REFF);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
    total += elapsed;
    cerr << "Time needed to compute all minimizers of Reference : ";
    cerr << (double)elapsed << "\n";
    cerr << SZ(retf) << " " << SZ(retr) << " " << SZ(for_index) << " " << SZ(rev_index) << "\n";
    ttp = fopen((char*)rr_str.c_str(), "r");
    fp = gzdopen(fileno(ttp), "r");
    seq = kseq_init(fp);
    int file_cnter = 0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    while (kseq_read(seq) >= 0)
    {
        string SS = string(seq->seq.s);
        string NN = string(seq->name.s);
        vector<string> strings;
        istringstream f(NN);
        string s;
        while (getline(f, s, '_'))
        {
            strings.emplace_back(s);
        }
        int offs = stoi(strings[5]);
        int strt = stoi(strings[1]);
//        for(auto t : strings)
//        {
//            cout<<t<<" ";
//        }
//        cout<<endl;
        // results.emplace_back(pool.enqueue([file_cnter, NN, SS]
        // {
        if (strings[2] != "unaligned")
        {
            find_kmers_of_read_in_refer((strings[4] == "F"), offs, strt, NN, SS);
            //return file_cnter;
            // }));
        }
        file_cnter++;
        //cerr<<file_cnter<<"\n";
    }
    kseq_destroy(seq);
    gzclose(fp);
    fclose(ttp);
    //for (auto &t : results) t.get();
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
    cerr << "Total time needed to do all!: " << elapsed << "\n";
    total += elapsed;

    return 0;
}
