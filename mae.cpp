#include <cstdio>
#include <cstdint>
#include <vector>
#include <cassert>
#include <bitset>
#include <vector>
#include <queue>
#include <map>
typedef uint8_t U8;
typedef uint16_t U16;
typedef uint32_t U32;
typedef uint64_t U64;

const std::string file_path = "./szdat2";
const int n = 100000000 / 4;
// const int n = 2869440;
// const int n = 10;
const int DEBUG = 1;

void dprintf(const char* format, ...) {
    if(!DEBUG) return;
    va_list args;
    va_start(args, format);
    vprintf(format, args); 
    va_end(args);
}

struct PredictTree {
    int depth;
    std::vector<int> table;
    PredictTree(int bitwidth = 16) {
        // need 2^(depth + 1) byte num
        depth = bitwidth;
        table.resize(1 << (depth + 1), 0);
    }
    int get_prob(int x, int len) {
        int now = 1;
        for(int i = len - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
        }
        int p1 = table[now * 2 + 1];
        int p2 = table[now * 2];
        // dprintf("get_prob: %d, %d, %d, %d\n", p1, p2, p1 + p2, 4095 * ((double) p1 / (p1 + p2)));
        assert(p1 >= 0 && p2 >= 0);
        if(p1 + p2 == 0) return 2048;
        return 4095 * ((double) p1 / (p1 + p2));
    }
    void add(int x) {
        int now = 1;
        for(int i = depth - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
            table[now]++;
        }
    }
    void del(int x) {
        int now = 1;
        for(int i = depth - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
            table[now]--;
        }
    }
};

struct SizedPredictTree {
    int depth;
    std::vector<int> table;
    std::deque<int> history;
    int maxsize;
    SizedPredictTree(int bitwidth = 16, int maxsize = 1024) {
        // need 2^(depth + 1) byte num
        depth = bitwidth;
        table.resize(1 << (depth + 1), 0);
        this->maxsize = maxsize;
    }
    int get_prob(int x, int len) {
        int now = 1;
        for(int i = len - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
        }
        int p1 = table[now * 2 + 1];
        int p2 = table[now * 2];
        // dprintf("get_prob: %d, %d, %d, %d\n", p1, p2, p1 + p2, 4095 * ((double) p1 / (p1 + p2)));
        assert(p1 >= 0 && p2 >= 0);
        if(p1 + p2 == 0) return 2048;
        return 4095 * ((double) p1 / (p1 + p2));
    }
    void add(int x) {
        history.push_back(x);
        int now = 1;
        for(int i = depth - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
            table[now]++;
        }
        if(history.size() > maxsize) {
            int y = history.front();
            history.pop_front();
            del(y);
        }
    }
    void del(int x) {
        int now = 1;
        for(int i = depth - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
            table[now]--;
        }
    }

};

struct Predictor {
    int pr, pr1, pr2;
    int w1, w2; // [0-10]
    SizedPredictTree* stree;
    PredictTree* tree;
    // use which tree to predict, current bitpos and last 32bit
    U32 bitpos, buffer;
    U32 bitwidth; 
    Predictor(int bitwidth = 16): pr(2048), bitpos(0), 
        buffer(0), bitwidth(bitwidth), 
        tree(new PredictTree(bitwidth)), stree(new SizedPredictTree(bitwidth, 1024)),
        w1(1), w2(1) {}

    bool update_bit(int y) {
        int ret = (pr > 2048 && y) || (pr < 2048 && !y);
        int r1 = (pr1 > 2048 && y) || (pr1 < 2048 && !y);
        int r2 = (pr2 > 2048 && y) || (pr2 < 2048 && !y);
        if(r1) {
            if(w1 > 3) w1 <<= 1;
            else w1++;
        } else {
            w1 = 1;
        }
        if(r2) {
            if(w2 > 3) w2 <<= 1;
            else w2++;
        } else {
            if(w2 > 1) w2--;
        }
        // dprintf("predict: %d, actual bit: %d\n", pr, y);
        bitpos++;
        buffer += buffer + y;
        if(bitpos == bitwidth) {
            tree->add(buffer); 
            stree->add(buffer);
            w1 = w2 = 1;
            bitpos = 0;       
            buffer = 0;
        }
        return ret;
    }
    int predict_bit() {
        int pr1 = stree->get_prob(buffer, bitpos);
        int pr2 = tree->get_prob(buffer, bitpos);
        pr = (pr1 * w1 + pr2 * w2) / (w1 + w2);
        // dprintf("predict: %d\n", pr);
        assert(pr >= 0 && pr < 4096);
        return pr;
    }
};

struct Encoder {
    U8* origin; int fb; 
    U32 x1, x2;    
    int bitwidth; 
    U32 correctcnt;
    Predictor* predictor;
    U32 x;      
    Encoder(int bitwidth = 16) {
        this->bitwidth = bitwidth;
        x1 = 0, x2 = 0xffffffff, x = 0;
        correctcnt = 0; fb = 0;
        predictor = new Predictor(bitwidth);
    }    

    int code(int y = 0) {
        // dprintf("codebit %d\n", y);
        int p = predictor->predict_bit();
        // dprintf("codebitx %d, p: %d\n", y, p);
        assert(p >= 0 && p < 4096);
        p += p < 2048;
        U32 xmid = x1 + ((x2 - x1) >> 12) * p + ((x2 - x1 & 0xfff) * p >> 12);
        assert(xmid >= x1 && xmid < x2);
        correctcnt += predictor->update_bit(y);
        y ? (x2 = xmid) : (x1 = xmid + 1);
        while (((x1 ^ x2) & 0xff000000) == 0) {
            origin[fb++] = x2 >> 24;
            x1 <<= 8;
            x2 = (x2 << 8) + 255;
        }
        return y;
    }
    int decode(int y = 0) {
        int p = predictor->predict_bit();
        assert(p >= 0 && p < 4096);
        p += p < 2048;
        U32 xmid = x1 + ((x2 - x1) >> 12) * p + ((x2 - x1 & 0xfff) * p >> 12);
        assert(xmid >= x1 && xmid < x2);
        correctcnt += predictor->update_bit(y);
        y = x <= xmid;
        y ? (x2 = xmid) : (x1 = xmid + 1);
        while (((x1 ^ x2) & 0xff000000) == 0) {
            x1 <<= 8;
            x2 = (x2 << 8) + 255;
            int cur = origin[fb++];
            x = (x << 8) + (cur & 255);
        }
        return y;
    }

    void compressbyte(U8 c) {
        for (int i = 7; i >= 0; --i)
            code((c >> i) & 1);
    }
    U8 decompressbyte() {
        int c = 0;
        for (int i = 0; i < 8; ++i)
            c += c + decode();
        return c;
    }
};

int bitwidth, base, basepos;
U32 successcnt = 0;

void compress(int *origin, size_t length, U8** res, int* res_len) {
    // dprintf("length: %ld\n", length);
    std::map<int, int> mp;
    bitwidth = 0;
    for(int i = 0; i < length; i++) {
        mp[origin[i]]++;
        // dprintf("origin[%d]: %d\n", i, origin[i]);
    }
    auto t = mp.end();
    t--;
    int maxx = t->first;
    int base = mp.begin()->first;
    unsigned int tmp = maxx - base;
    while(tmp) {
        bitwidth++;
        tmp >>= 1;
    }
    // if(base != 0) {
        for(int i = 0; i < length; i++) {
            origin[i] -= base;
        }
    // } else {
    //     int k = mp.begin()->first;
    //     int knum = mp.begin()->second;
    //     for(auto s: mp) {
    //         if(s.second > knum) {
    //             k = s.first;
    //             knum = s.second;
    //         }
    //     }
    //     // t == 0 -> t =  k - 1
    //     // t < k -> t = t - 1
    //     // t > k -> t = t
    //     for(int i = 0; i < length; i++) {
    //         if(!origin[i]) {
    //             origin[i] = k - 1;
    //         } else if(origin[i] < k) {
    //             origin[i]--;
    //         }
    //     }
    // }

    dprintf("bitwidth: %d, base: %d, maxx: %d\n", bitwidth, base, maxx);
    struct Encoder encoder(bitwidth);

    // actual compression process
    *res = (U8*)malloc((length * 2) * sizeof(U8));
    encoder.origin = *res;
    int last_correct = 0;
    for(int i = 0; i < length; i++) {
        // dprintf("code %d\n", origin[i]);
        for(int j = bitwidth - 1; j >= 0; j--) {
            encoder.code((origin[i] >> j) & 1);
        }
        if(i % 100000 == 0 && i) {
            dprintf("compress sucessful bit num: %d\n", encoder.correctcnt - last_correct);
            last_correct = encoder.correctcnt;
        }
    }
    *res_len = encoder.fb + 1;
    *res = (U8*)realloc(*res, *res_len);
    (*res)[encoder.fb] = encoder.x2 >> 24;
    successcnt  = encoder.correctcnt;
}

void decompress(U8 *origin, size_t length, int* res, int res_len) {
    struct Encoder encoder(bitwidth);
    encoder.origin = origin;
    for (int i = 0; i < 4; ++i) {
        encoder.x = (encoder.x << 8) + (origin[encoder.fb++] & 255);
    }
    for(int i = 0; i < res_len; i++) {
        res[i] = 0;
        for(int j = 0; j < bitwidth; j++) {
            res[i] += res[i] + encoder.decode();
        }
        res[i] += base;
    }
}

std::string debugU8ToString(U8 x) {
    std::string res = "";
    for(int i = 7; i >= 0; i--) {
        res += ((x >> i) & 1) + '0';
    }
    return res;
}

signed main() { 
    int* origin = (int*)malloc(n * sizeof(int));
    FILE *in = fopen(file_path.c_str(), "rb");
    for(int i = 1; i <= n; i++) {
        int s;
        fread(&s, sizeof(int), 1, in);
        origin[i - 1] = s;
    }
    fclose(in);
    U8* res;
    int res_len;
    compress(origin, n, &res, &res_len);
    dprintf("original file_size: %d," 
            "Compressed length: %d,"
            "total bit num: %lld,"
            "successful predict bit num: %lld\n", 
            n * 4, res_len, (U32)n * bitwidth, successcnt);
    return 0;
}
