#include <iostream>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cstdlib>
using namespace std;
using u64 = uint64_t;
using u32 = uint32_t;
using i64 = int64_t;
// Utilities 

static inline u64 nextpower(u64 x){     // compute next power of 2 for x
        if (x<=1) return 1;
        x--;
        x|=(x>>1);
        x|=(x>>2);
        x|=(x>>4);
        x|=(x>>8);
        x|=(x>>16);
        x|=(x>>32);
        x++;
        return x;
}
static inline int floorlog2(u64 x){      // compute floor of (log2(x))
        assert(x>0);
        return 63-__builtin_clzll(x);  
        // clz stands for "count leading zeros"     
}
static inline u64 highsqrt(u64 u){            //compute high bits of u (or quotient on division by root(u))
        int k=floorlog2(u);
        int half=(k+1)/2;
        return (u64)1<<half;
}
static inline u64 lowsqrt(u64 u){             //compute low bits of u (or remainder on division by root(u))
        int k=floorlog2(u);
        int half=k/2;
        return (u64)1<<half;
}
// Hashmap code to store clusters in the space-optimized vEB tree.
struct HMNode{
        u32 key;
        void* val;
        HMNode* next;
};
struct simplehashmap{
        size_t cap;
        HMNode** buckets;
        // Constructor
        simplehashmap(size_t initialcap=8){
                cap=initialcap;
                buckets=new HMNode*[cap];
                for (size_t i=0;i<cap;i++) buckets[i]=nullptr;
        }
        // Destructor
        ~simplehashmap(){
                for (size_t i=0;i<cap;i++){
                        HMNode* cur=buckets[i];
                        while(cur){
                                HMNode* nxt=cur->next;
                                delete cur;
                                cur=nxt;
                        }
                }
                delete[] buckets;
        }

        static inline u32 hash32(u32 x){
                x ^= x >> 16;
                x *= 0x7feb352d;
                x ^= x >> 15;
                x *= 0x846ca68b;
                x ^= x >> 16;
                return x;
        }
        void* get(u32 key) const{
                size_t idx=hash32(key)%cap;
                HMNode* cur=buckets[idx];
                while(cur){
                        if (cur->key==key) return cur->val;
                        cur=cur->next;
                }
                return nullptr;
        }
        bool contains(u32 key) const{
                return get(key)!=nullptr;
        }
        void put(u32 key,void* val){
                size_t idx=hash32(key)%cap;
                HMNode* cur=buckets[idx];
                while(cur){
                        if (cur->key==key){
                                cur->val=val;
                                return;
                        }
                        cur=cur->next;
                }
                HMNode* nn=new HMNode();
                nn->key=key;
                nn->val=val;
                nn->next=buckets[idx];
                buckets[idx]=nn;
        }
        bool erase(u32 key){
                size_t idx=hash32(key)%cap;
                HMNode* cur=buckets[idx];
                HMNode* prev=nullptr;
                while(cur){
                        if (cur->key==key){
                                if (prev) prev->next=cur->next;
                                else buckets[idx]=cur->next;
                                delete cur;
                                return true;
                        }
                        prev=cur;
                        cur=cur->next;
                }
                return false;
        }
};

// Optimized sparse veb
struct vebsparse {
        u64 u;
        i64 minv, maxv;
        vebsparse *summary;
        simplehashmap *clusters;

        //Optimization flags
        bool smallMode;         // bitset mode (u ≤ 64)
        bool mediumMode;        // byte array mode (64 < u ≤ 256)

        //Small mode storage
        u64 bitmask;            // up to 64 bits

        //Medium mode storage
        unsigned char *arr;     // size = u (<=256)

        vebsparse(u64 universe):
                u(nextpower(universe)),
                minv(-1),
                maxv(-1),
                summary(nullptr),
                clusters(nullptr),
                smallMode(false),
                mediumMode(false),
                bitmask(0),
                arr(nullptr)
        {
                if (u <= 64){
                        smallMode = true;
                        return;
                }
                if (u <= 256){
                        mediumMode = true;
                        arr = new unsigned char[u];
                        memset(arr, 0, u);
                        return;
                }
        }
        ~vebsparse(){
                if (smallMode) return;
                if (mediumMode){ delete[] arr; return; }

                if (summary) delete summary;
                if (clusters){
                        u64 ccount = highsqrt(u);
                        for (u64 i = 0; i < ccount; i++){
                                vebsparse *c = (vebsparse*)clusters->get((u32)i);
                                if (c) delete c;
                        }
                        delete clusters;
                }
        }

        inline bool empty() const { return minv == -1; }
        inline i64 minimum() const { return minv; }
        inline i64 maximum() const { return maxv; }

        inline u64 high(u64 x) const { u64 l = lowsqrt(u); return x / l; }
        inline u64 low(u64 x) const { u64 l = lowsqrt(u); return x & (l-1); }
        inline u64 index(u64 i, u64 j) const { u64 l = lowsqrt(u); return i*l + j; }

        // Small-mode functions (u ≤ 64)
        inline bool sm_member(u64 x) const{
                return (bitmask >> x) & 1ULL;
        }
        inline void sm_insert(u64 x){
                bitmask |= (1ULL << x);
        }
        inline void sm_erase(u64 x){
                bitmask &= ~(1ULL << x);
        }
        inline i64 sm_min() const{
                if (bitmask == 0) return -1;
                return __builtin_ctzll(bitmask);
        }
        inline i64 sm_max() const{
                if (bitmask == 0) return -1;
                return 63 - __builtin_clzll(bitmask);
        }
        inline i64 sm_succ(u64 x) const{
                u64 m = bitmask & (~((1ULL << (x+1)) - 1));
                if (m == 0) return -1;
                return __builtin_ctzll(m);
        }
        inline i64 sm_pred(u64 x) const{
                u64 mask = bitmask & ((1ULL << x) - 1);
                if (mask == 0) return -1;
                return 63 - __builtin_clzll(mask);
        }

        // Medium-mode functions (u ≤ 256)
        inline bool md_member(u64 x) const{
                return arr[x] != 0;
        }
        inline void md_insert(u64 x){
                arr[x] = 1;
        }
        inline void md_erase(u64 x){
                arr[x] = 0;
        }
        inline i64 md_min() const{
                for (u64 i = 0; i < u; i++)
                        if (arr[i]) return i;
                return -1;
        }
        inline i64 md_max() const{
                for (i64 i = (i64)u-1; i >= 0; i--)
                        if (arr[i]) return i;
                return -1;
        }
        inline i64 md_succ(u64 x) const{
                for (u64 i = x+1; i < u; i++)
                        if (arr[i]) return i;
                return -1;
        }
        inline i64 md_pred(u64 x) const{
                for (i64 i = (i64)x-1; i >= 0; i--)
                        if (arr[i]) return i;
                return -1;
        }
        //Main functions
        vebsparse* clusterof(u64 idx) const{
                if (!clusters) return nullptr;
                return (vebsparse*)clusters->get((u32)idx);
        }

        void putcluster(u64 idx, vebsparse *node){
                if (!clusters) clusters = new simplehashmap(16);
                clusters->put((u32)idx, node);
        }
        //INSERT
        void insertempty(u64 x){
                minv = maxv = (i64)x;
        }
        void insert(u64 x){
                //small-mode
                if (smallMode){
                        if (minv == -1){
                                sm_insert(x);
                                minv = maxv = x;
                                return;
                        }
                        sm_insert(x);
                        if ((i64)x < minv) minv = x;
                        if ((i64)x > maxv) maxv = x;
                        return;
                }
                //medium-mode
                if (mediumMode){
                        if (minv == -1){
                                md_insert(x);
                                minv = maxv = x;
                                return;
                        }
                        md_insert(x);
                        if ((i64)x < minv) minv = x;
                        if ((i64)x > maxv) maxv = x;
                        return;
                }
                //full vEB
                if (__builtin_expect(minv == -1, 0)){
                        insertempty(x);
                        return;
                }
                if ((i64)x < minv){
                        i64 temp = minv;
                        minv = x;
                        x = temp;
                }
                if (u > 2){
                        u64 h = high(x);
                        u64 l = low(x);
                        vebsparse *c = clusterof(h);

                        if (!c || c->empty()){
                                if (!summary) summary = new vebsparse(highsqrt(u));
                                if (!c){
                                        c = new vebsparse(lowsqrt(u));
                                        putcluster(h, c);
                                }
                                summary->insert(h);
                                c->insertempty(l);
                        } else {
                                c->insert(l);
                        }
                }
                if ((i64)x > maxv) maxv = x;
        }
        // ERASE
        void erase(u64 x){
                //small-mode
                if (smallMode){
                        if (minv == maxv){
                                bitmask = 0;
                                minv = maxv = -1;
                                return;
                        }
                        sm_erase(x);
                        minv = sm_min();
                        maxv = sm_max();
                        return;
                }
                //medium-mode
                if (mediumMode){
                        if (minv == maxv){
                                memset(arr, 0, u);
                                minv = maxv = -1;
                                return;
                        }
                        md_erase(x);
                        minv = md_min();
                        maxv = md_max();
                        return;
                }
                //full vEB 
                if (minv == maxv){
                        minv = maxv = -1;
                        return;
                }
                if (u == 2){
                        minv = (x == 0 ? 1 : 0);
                        maxv = minv;
                        return;
                }
                if ((i64)x == minv){
                        u64 fc = summary->minimum();
                        vebsparse *c = clusterof(fc);
                        u64 newlow = c->minimum();
                        u64 newx = index(fc, newlow);
                        minv = newx;
                        x = newx;
                }
                u64 h = high(x), l = low(x);
                vebsparse *c = clusterof(h);
                c->erase(l);
                if (c->empty()){
                        clusters->erase((u32)h);
                        delete c;
                        summary->erase(h);
                        if (summary->empty()){
                                delete summary;
                                summary = nullptr;
                        }
                }
                if ((i64)x == maxv){
                        if (!summary){
                                maxv = minv;
                        } else {
                                u64 mc = summary->maximum();
                                vebsparse *c2 = clusterof(mc);
                                maxv = index(mc, c2->maximum());
                        }
                }
        }
        //MEMBER 
        bool member(u64 x) const{
                if (smallMode) return sm_member(x);
                if (mediumMode) return md_member(x);
                if (minv == (i64)x || maxv == (i64)x) return true;
                if (u == 2) return false;
                u64 h = high(x);
                vebsparse *c = clusterof(h);
                if (!c) return false;
                return c->member(low(x));
        }
        //SUCCESSOR 
        i64 successor(u64 x) const{
                if (smallMode) return sm_succ(x);
                if (mediumMode) return md_succ(x);
                if (u == 2){
                        if (x == 0 && maxv == 1) return 1;
                        return -1;
                }
                if (minv != -1 && x < (u64)minv) return minv;
                u64 h = high(x), l = low(x);
                vebsparse *c = clusterof(h);
                if (c && c->maximum() != -1 && l < (u64)c->maximum()){
                        i64 off = c->successor(l);
                        return index(h, off);
                } else {
                        if (!summary) return -1;
                        i64 sc = summary->successor(h);
                        if (sc == -1) return -1;
                        vebsparse *c2 = clusterof(sc);
                        return index(sc, c2->minimum());
                }
        }
        //PREDECESSOR 
        i64 predecessor(u64 x) const{
                if (smallMode) return sm_pred(x);
                if (mediumMode) return md_pred(x);
                if (u == 2){
                        if (x == 1 && minv == 0) return 0;
                        return -1;
                }
                if (maxv != -1 && x > (u64)maxv) return maxv;
                u64 h = high(x), l = low(x);
                vebsparse *c = clusterof(h);
                if (c && c->minimum() != -1 && l > (u64)c->minimum()){
                        i64 off = c->predecessor(l);
                        return index(h, off);
                } else {
                        if (!summary){
                                if (minv != -1 && x > (u64)minv) return minv;
                                return -1;
                        }
                        i64 pc = summary->predecessor(h);
                        if (pc == -1){
                                if (minv != -1 && x > (u64)minv) return minv;
                                return -1;
                        }
                        vebsparse *c2 = clusterof(pc);
                        return index(pc, c2->maximum());
                }
        }
};
static void test_sparse(){
        cout << "Enter universe size :";
        int u;
        cin >> u;
        string operations[]={"Exit","Insert","Delete","Successor","Predecessor","Maximum","Minimum","Member"};
        vebsparse v(u);
        cout << "Menu: " << endl;
        for (int i=0;i<8;i++){
                cout << i << '=' << operations[i] << ' ';
        }
        cout << endl;
        int ch;
        cout << "Enter your choice:";
        cin >> ch;
        while(ch!=0){
                if (ch==1){
                        int x;
                        cout << "Enter the number to insert:";
                        cin >> x;
                        v.insert(x);
                        cout << "Successfully inserted" << endl;
                }
                else if (ch==2){
                        int x;
                        cout << "Enter the number to delete:";
                        cin >> x;
                        if (v.member(x)){
                                v.erase(x);
                                cout << "Successfully deleted" << endl;
                        }
                        else cout << "No such element exists" << endl;
                }
                else if (ch==3){
                        int x;
                        cout << "Successor of:";
                        cin >> x;
                        int y=v.successor(x);
                        if (y==-1) cout << "No successor" << endl;
                        else cout << y << endl;
                }
                else if (ch==4){
                        int x;
                        cout << "Predecessor of:";
                        cin >> x;
                        int y=v.predecessor(x);
                        if (y==-1) cout << "No predecessor" << endl;
                        else cout << y << endl;
                }
                else if (ch==5){
                        int y=v.maximum();
                        if (y==-1) cout << "No element present" << endl;
                        else cout << y << endl;  
                }
                else if (ch==6){
                        int y=v.minimum();
                        if (y==-1) cout << "No element present" << endl;
                        else cout << y << endl;  
                }
                else if (ch==7){
                        int y;
                        cout << "Enter the number to search:";
                        cin >> y;
                        if (v.member(y)) cout << "Exists in tree" << endl;
                        else cout << "Does not exist in tree" << endl;
                }
                else cout << "Invalid choice" << endl;
                cout << "Enter your choice:";
                cin >> ch;  
        }
} 
int main(){
        test_sparse();
        return 0;
}