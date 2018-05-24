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

    void repeatPixel(int x, int y) {
        if (x==0) return;
        auto i = y*W + x;
        colors[i] = colors[i-1];
    }

    void clear() {
        colors[] = 0;
        depth[] = 1000000.0;
    }

    void fillRows(int y0, int h, uint clr) {
        colors[y0*W .. (y0+h)*W] = clr;
    }
}

version(DirectorsCut) {
// pragma(msg, "img: DirectorsCut");
import dimage.stream, std.stdio;

void saveBuf(string fname, ColorDrawBuf cdbuf) {
    import dimage.png, dimage.image, std.parallelism;
    const W = cdbuf.width, H = cdbuf.height;
    // static SuperImage simg;
    // if (simg is null)
    auto simg = new SuperImage(W, H, 3, 8);
    foreach(y; 0..H) {
        uint* ptr = cdbuf.scanLine(y);
        foreach(x; 0..W)
            simg[x,y] = ptr[x];
    }
    taskPool.put( task({
        OutputStream output = new SimpleOutStream(fname);
        savePNG(simg, output);
        output.close();
    }));
}

class SimpleOutStream : OutputStream {
    File f;

    this(string fname) { f = File(fname, "wb"); }

    override ulong getPosition() @property {
        return f.tell();
    }
    override bool setPosition(ulong pos) {  f.seek(pos); return true; }
    override ulong size() { return f.size; }
    override void close() { return f.close(); }
    override bool seekable() { return true; }
    override void flush() { f.flush(); }
    override bool writeable() { return true; }
    override ulong writeBytes(const(void*) buffer, ulong count) {
        f.rawWrite(buffer[0..count]); return count;
    }
}

}