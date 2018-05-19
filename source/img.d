import dlangui;

class ImageZ {
    int W,H;
    uint[] colors;
    double[] depth;

    this(int w, int h) {
        W = w; H = h;
        colors.length = W*H;
        depth.length = W*H;
        depth[] = 1000000.0;
    }

    void copyFrom(ImageZ img) {
        assert(W==img.W && H==img.H);
        colors[] = img.colors[];
        depth[] = img.depth[];
    }

    void putPixelZ(int x, int y, uint clr, double z) {
        if (!(x >= 0 && x < W)) return;
        if (!(y >=0 && y < H)) return;
        auto i = y*W + x;
        if (z < depth[i]) {
            depth[i] = z;
            colors[i] = clr;
        }
    }

    void putPixel(int x, int y, uint clr) {
        if (!(x >= 0 && x < W)) return;
        if (!(y >=0 && y < H)) return;
        auto i = y*W + x;
        colors[i] = clr;
    }

    void clear() {
        colors[] = 0;
        depth[] = 1000000.0;
    }

    void fillRows(int y0, int h, uint clr) {
        colors[y0*W .. (y0+h)*W] = clr;
    }
}
