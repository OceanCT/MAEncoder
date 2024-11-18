#include <cstdio>
#include <cstdint>
#include <vector>
#include <cassert>
#include <bitset>
#include <vector>
#include <queue>
#include <map>
#include <unistd.h>
#include <cmath>
typedef uint8_t U8;
typedef uint16_t U16;
typedef uint32_t U32;
typedef uint64_t U64;
const std::string file_path = "./szdat2";
// const std::string file_path = "./data_pre/CESM-ATM/szdats/CLOUD_1_26_1800_3600.f32.sz4.szdat";
// const int n = 100000000 / 4;
const int n = 100000 * 10;
// const int n = 673920000 / 4;
// const int n = 2869440;
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
    int curdepth;
    int pos; // [0, depth]
    std::vector<int> table;
    PredictTree(int bitwidth = 16) {
        // need 2^(depth + 1) byte num
        depth = bitwidth;
        table.resize(1 << (depth + 1), 0);
        pos = 1;
        curdepth = 0;
    }
    int get_currentprob() {
        int p1 = table[pos * 2 + 1];
        int p2 = table[pos * 2];
        assert(p1 >= 0 && p2 >= 0);
        if(p1 + p2 == 0) return 2048;
        if(p1 <= 50 && p2 <= 50) p1++, p2++;
        double s = (double) (p1 - p2) / (p1 + p2);
        int pr1 = 0, pr2 = 0, pr3 = 0, pr05 = 0;
        pr1 = 2047 + 2048 * s;
        if(p1 > p2) pr2 = 2047 + 2048 * s * s;
        else pr2 = 2047 - 2048 * s * s;
        pr3 = 2047 + 2048 * s * s * s;
        if(p1 > p2) pr05 = 2047 + 2048 * sqrt(s);
        else pr05 = 2047 - 2048 * sqrt(-s);
        int pr = pr1;
        if(pr >= 4096) pr = 4095;
        if(pr < 0) pr = 0;
        // printf("p1: %d, p2: %d, pr: %d\n", p1, p2, pr);
        // return pr;
        return pr;
        // return 4095 * ((double) p1 / (p1 + p2));
    }
    void updatebit(int y) {
        pos = pos * 2 + y;
        table[pos]++;
        curdepth++;
        if(curdepth == depth) {
            curdepth = 0;
            pos = 1;
        }
    }
};

struct SizedPredictTree {
    int depth;
    std::vector<int> table;
    std::deque<int> history;
    int curdepth;
    int pos; // [0, depth]
    int maxsize;
    SizedPredictTree(int bitwidth = 16, int maxsize = 1024) {
        depth = bitwidth;
        table.resize(1 << (depth + 1), 0);
        this->maxsize = maxsize;
        pos = 1;
        curdepth = 0;
    }
    int get_currentprob() {
        int p1 = table[pos * 2 + 1];
        int p2 = table[pos * 2];
        // printf("p1: %d, p2: %d\n", p1, p2);
        assert(p1 >= 0 && p2 >= 0);
        if(p1 + p2 == 0) return 2048;
        if(p1 <= 50 && p2 <= 50) p1++, p2++;
        if(!p1 || !p2) p1++, p2++;
        return 4095 * ((double) p1 / (p1 + p2));
    }
    void add(int x) {
        history.push_back(x);
        if(history.size() > maxsize) {
            int y = history.front();
            history.pop_front();
            del(y);
        }
    }
    void updatebit(int x) {
        pos = pos * 2 + x;
        table[pos]++;
        curdepth++;
        if(curdepth == depth) {
            curdepth = 0;
            pos = 1;
        }
    }
    void del(int x) {
        int now = 1;
        for(int i = depth - 1; i >= 0; i--) {
            now = now * 2 + ((x >> i) & 1);
            table[now]--;
            assert(table[now] >= 0);
        }
    }

};


int printcnt = 0;
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
    void init() {
        stree->depth = bitwidth;
        tree->depth = bitwidth;
    }
    bool update_bit(int y) {
        stree->updatebit(y);
        tree->updatebit(y);
        int ret = (pr > 2048 && y) || (pr < 2048 && !y);
        bool r1 = (pr1 > 2048 && y) || (pr1 < 2048 && !y);
        bool r2 = (pr2 > 2048 && y) || (pr2 < 2048 && !y);
        // printf("bitpos: %d, act: %d, pr: %d, w1: %d, w2: %d\n", bitpos, y, pr, w1, w2);
        if(r1) {
            if(w1 > 4) w1 <<= 1;
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
        bitpos++;
        buffer += buffer + y;
        if(bitpos == bitwidth) {
            stree->add(buffer);
            // printf("\n");
            // for(int sss = bitwidth - 1; sss >= 0; sss--) {
            //     printf("%d,", (buffer >> sss) & 1);
            // }
            w1 = w2 = 1;
            bitpos = 0;       
            buffer = 0;
        }
        return ret;
    }
    int predict_bit() { 
        pr1 = stree->get_currentprob();
        pr2 = tree->get_currentprob();
        pr = (pr1 * w1 + pr2 * w2) / (w1 + w2);
        assert(pr >= 0 && pr < 4096);
        return pr;
    }
};

int cnt = 0;
int bitwidth = 0;

struct Encoder {
    U8* origin; int fb; 
    U32 x1, x2, x;    
    int bitwidth; 
    U32 correctcnt;
    struct Predictor predictor;  
    Encoder(int bitwidth = 32) {
        this->bitwidth = bitwidth;
        x1 = 0, x2 = 0xffffffff, x = 0;
        correctcnt = 0; fb = 0;
        predictor.bitwidth = bitwidth;
        predictor.init();
    }    

    int code(int y = 0) {
        int p = predictor.predict_bit();
        assert(p >= 0 && p < 4096);
        p += p < 2048;
        U32 xmid = x1 + ((x2 - x1) >> 12) * p + ((x2 - x1 & 0xfff) * p >> 12);
        assert(xmid >= x1 && xmid < x2);
        y ? (x2 = xmid) : (x1 = xmid + 1);
        correctcnt += predictor.update_bit(y);
        while (((x1 ^ x2) & 0xff000000) == 0) {
            origin[fb++] = x2 >> 24;
            x1 <<= 8;
            x2 = (x2 << 8) + 255;
        }
        // printf("file size: %d\n", fb);
        return y;
    }
    int decode(int y = 0) {
        int p = predictor.predict_bit();
        assert(p >= 0 && p < 4096);
        p += p < 2048;
        U32 xmid = x1 + ((x2 - x1) >> 12) * p + ((x2 - x1 & 0xfff) * p >> 12);
        assert(xmid >= x1 && xmid < x2);
        y = x <= xmid;
        y ? (x2 = xmid) : (x1 = xmid + 1);
        correctcnt += predictor.update_bit(y);
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

U32 successcnt = 0;
void compress(int *origin, size_t length, U8** res, int* res_len) {
    bitwidth = 0;
    int maxx = origin[0];
    for(int i = 0; i < length; i++) {
        if(origin[i] > maxx) {
            maxx = origin[i];
        }
    }
    while(maxx) {
        bitwidth++;
        maxx >>= 1;
    }
    int reflection = (maxx + 1) / 2;
    // maxx = 5 reflct = 4
    // 0 1 2 3 4 5
    // 0 3 2 4 1 5
    // 0 reflect reflect - 1  reflect + 1

    dprintf("bitwidth: %d\n", bitwidth);
    // actual compression process
    struct Encoder encoder(bitwidth);
    *res = (U8*)malloc((length * 2) * sizeof(U8));
    (*res)[encoder.fb++] = bitwidth;
    encoder.origin = *res;
    int lastsuccess = 0;
    int lastsize = 0;
    for(int i = 0; i < length; i++) {
        if(i > 0 && i % 100000 == 0) {
            printf("successful cnt: %d\n", encoder.correctcnt - lastsuccess);
            printf("size: %d\n", encoder.fb - lastsize);
            lastsuccess = encoder.correctcnt;
            lastsize = encoder.fb;
        }
        if(origin[i] >= reflection) {
            origin[i] = (origin[i] - reflection) * 2 + 1;
        } else {
            origin[i] = (reflection - origin[i]) * 2;
        }
        for(int j = bitwidth - 1; j >= 0; j--) {
            encoder.code((origin[i] >> j) & 1);
        }
    }
    *res_len = encoder.fb + 1;
    *res = (U8*)realloc(*res, *res_len);
    (*res)[encoder.fb] = encoder.x2 >> 24;
    successcnt  = encoder.correctcnt;
}

void decompress(U8 *origin, size_t length, int* res, int res_len) {
    bitwidth = origin[0];
    struct Encoder encoder(bitwidth);
    dprintf("bitwidth: %d\n", bitwidth);
    encoder.origin = origin;
    encoder.fb = 1;
    for (int i = 0; i < 4; ++i) {
        encoder.x = (encoder.x << 8) + (origin[encoder.fb++] & 255);
    }
    for(int i = 0; i < res_len; i++) {
        res[i] = 0;
        for(int j = 0; j < bitwidth; j++) {
            res[i] += res[i] + encoder.decode();
        }
    }
}

int Mode = 1; // 1 for compress, 0 for decompress

signed main() { 
    int* origin = (int*)malloc(n * sizeof(int));
    FILE *in = fopen(file_path.c_str(), "rb");
    for(int i = 1; i <= n; i++) {
        int s;
        fread(&s, sizeof(int), 1, in);
        origin[i - 1] = s;
    }
    fclose(in);
    if(Mode) {
      U8* res;
      int res_len;
      time_t start = time(0);
      compress(origin, n, &res, &res_len);
      double seconds_since_start = difftime(time(0), start);
      printf("original file_size: %d," 
              "Compressed length: %d,"
              "total bit num: %u,"
              "successful predict bit num: %u,"
              "bitwidth: %d,"
              "time cost: %f\n", 
              n * 4, res_len, (U32)n * bitwidth, 
              successcnt, bitwidth,
              seconds_since_start);
      FILE* ft = fopen("mae.tmp", "wb");
      for(int i = 0; i < res_len; i++) {
          putc(res[i], ft);
      }
    } else {
      FILE* ft = fopen("mae.tmp", "rb");
      U8* res;
      int res_len;
      fseek(ft, 0, SEEK_END);
      res_len = ftell(ft);
      fseek(ft, 0, SEEK_SET);
      res = (U8*)malloc(res_len * sizeof(U8));
      for(int i = 0; i < res_len; i++) {
          res[i] = getc(ft);
      }
      int* origin1 = (int*)malloc(n * sizeof(int));
      decompress(res, res_len, origin1, n);
      for(int i = 0; i < n; i++) {
            // printf("origin[%d]: %d, origin1[%d]: %d\n", i, origin[i], i, origin1[i]);
            if(origin[i] != origin1[i]) {
                printf("origin[%d]: %d, origin1[%d]: %d\n", i, origin[i], i, origin1[i]);
                break;
            }
      }
    }
    return 0;
}
