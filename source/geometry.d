import std.math, std.array, std.format, std.stdio;
import img, symbolic, dlangui.graphics.drawbuf;
import std.algorithm : canFind, map, sum, min;
import std.range : iota;

__gshared uint[] wallTexture;
__gshared ColorDrawBuf worldMap;

struct UV { double u,v; }

// f'(t,y) = drv(y)
T evolveRK(alias drv, T)(T state, double dt) {
    auto a = drv(state);
    auto b = drv(state + a * (dt/2));
    auto c = drv(state + b * (dt/2));
    auto d = drv(state + c*dt);
    return state + a*(dt/6) + b*(dt/3) + c*(dt/3) + d*(dt/6);
}

struct Vec(T, int N) {
    T[N] data;

    Vec!(T,N) opBinary(string op)(double k) if (op=="*") {
        Vec!(T,N) v;
        v.data[] = data[] * k;
        return v;
    }

    Vec!(T,N) opBinary(string op)(Vec!(T,N) b) if (op=="+") {
        Vec!(T,N) v;
        v.data[] = data[] + b.data[];
        return v;
    }
}

alias Dbl4 = Vec!(double, 4);
enum wallColor = 0xF08000;
__gshared double globalR = 1000.0, globalDistance = 3.0;

string[] genCode(alias surfaceEquation)() {
    Expr[] X = surfaceEquation();
    string embed = format("x = %s; y = %s; z = %s;", X[0].code, X[1].code, X[2].code).addTrigNeeded;

    auto Xu = X.diff("u"); // basis vectors
    auto Xv = X.diff("v");
    auto Xus = format("Xu: %s", Xu);
    auto Xvs = format("Xv: %s", Xv);
    Expr[2][2] g, ginv; // metric tensor and its inverse
    g[0][0] = dot(Xu, Xu);
    g[0][1] = g[1][0] = dot(Xu, Xv);
    g[1][1] = dot(Xv, Xv);
    auto det_g = add(mul(g[0][0], g[1][1]), neg(mul(g[0][1], g[1][0])));
    ginv[0][0] = div(g[1][1], det_g);
    ginv[0][1] = ginv[1][0] = div(neg(g[0][1]), det_g);
    ginv[1][1] = div(g[0][0], det_g);

    Expr[2][2][2] C; // Christoffel symbols
    string[] dx = ["u","v"];
    auto two = new Const("2");
    foreach(k; 0..2)
        foreach(i; 0..2)
            foreach(j; 0..2) {
                //Ckij = 1/2 * ginv_kl * (d/dx_j g_li  + d/dx_i g_lj - d/dx_l g_ij )
                Expr[2] tmp;
                foreach(l; 0..2) {
                    auto J = g[l][i].diff(dx[j]);
                    auto I = g[l][j].diff(dx[i]);
                    auto L = neg(g[i][j].diff(dx[l]));
                    tmp[l] = div( mul(ginv[k][l], add(J, I, L)), two);
                }
                C[k][i][j] = add(tmp[0], tmp[1]);
            }
    // res.data[1] = -C[0][0][0]*du*du - 2 * C[0][0][1]*du*dv - C[0][1][1]*dv*dv; //u''
    // res.data[3] = -C[1][0][0]*du*du - 2 * C[1][0][1]*du*dv - C[1][1][1]*dv*dv; //v''
    string u_accel = format("res.data[1] = -%s*du*du - 2 * %s*du*dv - %s*dv*dv;\n",
                            C[0][0][0].code, C[0][0][1].code, C[0][1][1].code).txtSimp;
    string v_accel = format("res.data[3] = -%s*du*du - 2 * %s*du*dv - %s*dv*dv;\n",
                            C[1][0][0].code, C[1][0][1].code, C[1][1][1].code).txtSimp;
    string accel = addTrigNeeded(u_accel ~ v_accel, false);

    Expr du = new Var("du"), dv = new Var("dv");
    Expr mag = add(mul(g[0][0], du, du), mul(two, g[0][1], du, dv), mul(g[1][1], dv,dv));
    string vlen = format("return sqrt(%s);", mag.code).txtSimp.addTrigNeeded;
    return [embed, accel,  vlen, g[0][0].code, g[0][1].code, g[1][1].code, Xus, Xvs];
    //      0      1        2         3           4            5
}

string addTrigNeeded(string e, bool includeCosV=true) {
    return
        (e.canFind("sin_v") ? "const double sin_v = sin(v);\n" : "") ~
        (e.canFind("cos_v") && includeCosV ? "const double cos_v = cos(v);\n" : "") ~
        (e.canFind("sin_u") ? "const double sin_u = sin(u);\n" : "") ~
        (e.canFind("cos_u") ? "const double cos_u = cos(u);\n" : "") ~
        e;
}

uint earthColor(double u, double v) {
    int iu = cast(int)(u/PI*360), iv = cast(int)(v/PI*360);
    if (iu >= 720) {
        do iu -= 720; while (iu >= 720);
    } else if (iu < 0) {
        do iu += 720; while (iu < 0);
    }
    if (iv >= 360) {
        do iv -= 360; while (iv >= 360);
    } else if (iv < -360) {
        do iv += 360; while (iv < -360);
    }
    if (iv >= 180) {
        iv = 360 - iv;
    } else if (iv <= -180) {
        iv = -359 - iv;
    }
    iv = 180 - iv;
    return worldMap.scanLine(iv)[iu];
}

uint synthColor(double u, double v) {
    int iu = cast(int)(u/PI*180), iv = cast(int)(v/PI*180);
    if (iu >= 360) {
        do iu -= 360; while (iu >= 360);
    } else if (iu < 0) {
        do iu += 360; while (iu < 0);
    }
    if (iv >= 180) {
        do iv -= 360; while (iv >= 180);
    } else if (iv < -180) {
        do iv += 360; while (iv < -180);
    }

    if (iv >= 90) {
        iv = 180 - iv;
    } else if (iv < -90) {
        iv = -180 - iv;
    }
    // iu: 0..360,  iv: -90..90

    uint clr;
    int mu = iu % 20, mv = (iv + 180) % 20;
    if (mu >= 10 && mv >= 10 && mu < 16 && mv < 16)
        clr = wallColor; //squares
    else {
        int c = ((iu ^ iv) & 31)*4 + 64;
        uint rmask = iu >= 180 ? 255 : 127;
        uint gmask = ((iu / 90) & 1) ? 255 : 127;
        uint bmask = iv >= 0 ? 255 : 127;
        clr = ((c & rmask) << 16) + ((c & gmask)<<8) + (c & bmask);
    }
    return clr;
}



class Surface(alias surfaceEquation, bool poleSingularity) {
    enum codes = genCode!surfaceEquation();
    pragma(msg, codes[6]);
    pragma(msg, codes[7]);

    static void embedIn3D(double u, double v, ref double x, ref double y, ref double z) {
        pragma(msg, codes[0]);
        const double R = globalR;
        mixin(codes[0]);
    }

    static Dbl4 geodesicStep(Dbl4 y) {
        //y : u, u', v, v'
        Dbl4 res;
        auto du = res.data[0] = y.data[1]; // u' = u'
        auto dv = res.data[2] = y.data[3]; // v' = v'
        const double R = globalR;
        const double v = y.data[2], u = y.data[0];
        const double cos_v = cos(v);
        static if (poleSingularity)
            if (abs(cos_v) < 0.000001) { // handle poles, poorly
                res.data[1] = 0;
                res.data[3] = 0;
                return res;
            }
        pragma(msg, codes[1]);
        mixin(codes[1]);
        return res;
    }

    static double vlen(double u, double v, double du, double dv) {
        pragma(msg, codes[2]);
        const double R = globalR;
        mixin(codes[2]);
    }

    static void courseVector(double u, double v, double angle, ref double du, ref double dv) {
        const double R = globalR;
        double cos_v = cos(v);
        const double cos_u = cos(u), sin_u = sin(u), sin_v = sin(v);
        static if (poleSingularity)
            if (abs(cos_v) < 0.000001) cos_v = cos_v >= 0 ? 0.000001 : 0.000001;
        enum u_code = format("du = sin(angle) / sqrt(%s);", codes[3]).txtSimp;
        enum v_code = format("dv = cos(angle) / sqrt(%s);", codes[5]).txtSimp;
        pragma(msg, u_code);
        pragma(msg, v_code);
        mixin(u_code);
        mixin(v_code);
        //du = sin(angle)/(cos_v*R); dv = cos(angle)/R;
    }

    static double courseAngle(double u, double v, double du, double dv) { // in radians
        // up (0, 1/sqrt(g22=R^2)) = (0, 1/R) for sphere
        // up*V = cos(a) * |V| = g12 * up_v * du + g22 * up_v * dv = R * dv
        // cos(a) = R*dv / |V|
        const double R = globalR;
        const double cos_v = cos(v), sin_v = sin(v);
        const double sin_u = sin(u), cos_u = cos(u);
        enum up_v_code = format("double up_v = 1/sqrt(%s);", codes[5]).txtSimp;
        enum product_code = format("double product = %s * up_v * du + %s * up_v * dv;", codes[4], codes[5]).txtSimp;
        pragma(msg, up_v_code);
        pragma(msg, product_code);
        mixin(up_v_code);
        mixin(product_code);
        double a = acos(product / vlen(u, v, du, dv));
        if (du < 0) a = 2*PI - a;
        return a;
    }
}// Surface

Expr[] sphereEq() {
    auto R = new Var("R");
    return [mul(R, mul(new Cos("u"), new Cos("v"))),
            mul(R, new Sin("v")),
            mul(R, mul(new Sin("u"), new Cos("v")))  ];
}

Expr[] planeEq() {
    auto R = new Var("R");
    return [mul(new Var("u"), R), zero, mul(new Var("v"),R)];
}

Expr[] ellipsoidEq() {
    auto R = new Var("R");
    return [mul(R, mul(new Cos("u"), new Cos("v"))),
            mul(div(R, new Const("2")), new Sin("v")),
            mul(R, mul(new Sin("u"), new Cos("v")))  ];
}

Expr[] torusEq() {
    auto R = new Var("R"), two = new Const("2");
    return [mul(R, add(two, new Cos("v")), new Cos("u")),
            mul(R, add(two, new Cos("v")), new Sin("u")),
            mul(neg(R), new Sin("v")) ];
}

double ballPointColor(double z) { return -z/globalR; }

struct Sphere {
    static string name = "Ball";
    alias equation = sphereEq;
    static double distance = 3.0;
    alias color = synthColor;
    alias pointColor = ballPointColor;
    static bool singularity = true;
}

struct Earth {
    static string name = "Earth";
    alias equation = sphereEq;
    static double distance = 3.0;
    alias color = earthColor;
    alias pointColor = ballPointColor;
    static bool singularity = true;
}

struct Plane {
    static string name = "Plane";
    alias equation = planeEq;
    static double distance = 3.0;
    alias color = synthColor;
}

struct FatBall {
    static string name = "Fat ball";
    alias equation = ellipsoidEq;
    static double distance = 3.0;
    alias color = synthColor;
    alias pointColor = ballPointColor;
    static bool singularity = true;
}

struct Donut {
    static string name = "Donut";
    alias equation = torusEq;
    static double distance = 9.0;
    static const int NV = 700;
    static double pickV(int iv, int nv) {
        return (iv - nv/2) * PI*2 / nv;
    }

    static uint color(double u, double v) { // u: 0 .. 2Pi,   v: -Pi .. Pi
        int iu = cast(int)(u/PI*180), iv = cast(int)(v/PI*180);
        if (iu >= 360) {
            do iu -= 360; while (iu >= 360);
        } else if (iu < 0) {
            do iu += 360; while (iu < 0);
        }
        if (iv >= 360) {
            do iv -= 360; while (iv >= 360);
        } else if (iv < 0) {
            do iv += 360; while (iv < 0);
        }

        // iu: 0..360,  iv: 0..360
        uint clr;
        int mu = iu % 20, mv = (iv + 180) % 20;
        if (mu >= 10 && mv >= 10 && mu < 16 && mv < 16)
            clr = wallColor; //squares
        else {
            int c = ((iu ^ iv) & 31)*4 + 64;
            uint rmask = iu >= 180 ? 255 : 127;
            uint gmask = ((iu / 90) & 1) ? 255 : 127;
            uint bmask = iv >= 180 ? 255 : 127;
            clr = ((c & rmask) << 16) + ((c & gmask)<<8) + (c & bmask);
        }
        return clr;
    }

    static double pointColor(double z) { return -z/(globalR*3.0); }
}//Donut

struct BlackHole {
    static string name = "Black hole";
    static Expr[] equation() {
        auto R = new Var("R");
        auto r = mul(new Var("v"), R);
        return [mul(r, new Cos("u")),
                mul(R, new Fv()),
                mul(r, neg(new Sin("u"))) ];

    }
    static const int NV = 700;
    static double pickV(int iv, int nv) {
        return 1.0 + iv * 8.0 / nv;  // 1..9
    }
    static double distance = 12.0;
    static uint color(double u, double v) { // u: 0 .. 2Pi,   v: 1..9
        const double R = globalR;
        int x = cast(int) (v*cos(u)*R)/10;
        int z = cast(int) (-v*sin(u)*R)/10;
        int c = ((x ^ z) & 31)*4 + 64;
        uint rmask = x > 0 ? 255 : 127;
        uint gmask = z > 0 ? 255 : 127;
        uint bmask = ((x/64 + z/64)&1) ? 255 : 127;
        return ((c & rmask) << 16) + ((c & gmask)<<8) + (c & bmask);
    }
    static double pointColor(double z) { return 1.0; }
    static void ensureCorrectPos(Params ps) {
        if (ps.pos.v < 4.0) ps.pos.v = 4.0;
    }
}

// In polar coordinates where u is angle and v is radius and height is F(v)
// metric is ds^2 = (1 + dF2(v)) * dv^2 + v^2 * du^2
// while Schwarzschild metric, just the space part, is
// ds^2 = (1 + R/(r-R)) * dr^2 + r^2 * du^2
// where R is Schwarzschild radius
// so dF2(r) = R/(r-R)
// dF(r) = sqrt(R/(r-R))
// F(r) = 2*sqrt(R*(r-R))
// dF(r) * ddF = -R / (2*(R-r)**2)
// C_111 = dF * ddF / (1 + dF*dF) = R / (2*r*(R-r))
// here R = 1
double F(double r) {  return 2.0 * sqrt(r - 1.0); }
double dF2(double r) { return 1.0 / (r-1.0); }
double GC111(double r) { return 0.5 / (r*(1.0-r)); }

struct Wormhole {
    static string name = "Wormhole";
    static Expr[] equation() {
        auto R = new Var("R");
        auto v = new Var("v");
        auto r = mul(add(mul(v, v), new Const("1")),  R);
        return [mul(r, new Cos("u")),
                mul(r, new Sin("u")),
                mul(v, R)];
    }
    static const int NV = 1200;
    static double pickV(int iv, int nv) {
        return (iv - nv/2) * 10.0 / nv;  // -5..5
    }
    static double distance = 12.0;
    static uint color(double u, double v) { // u: 0 .. 2Pi,   v: -5..5
        const double R = globalR;
        auto r = (1 + v*v)*R;
        int x = cast(int) (r*cos(u))/10;
        int z = cast(int) (r*sin(u))/10;
        int c = ((x ^ z) & 31)*4 + 64;
        uint rmask = x > 0 ? 255 : 127;
        uint gmask = z > 0 ? 255 : 127;
        uint bmask = v > 0 ? 255 : 127;
        return ((c & rmask) << 16) + ((c & gmask)<<8) + (c & bmask);
    }
    static double pointColor(double z) { return 1.0; }
}

class Params {
    double heading, dt, range;
    UV pos;
    bool dyndt, walls;
    double rotAlpha, rotBeta; // rotation angles in degrees
}

double[3][3] rotMatrix(Params ps) {
    double alpha = ps.rotAlpha * PI / 180, beta = ps.rotBeta * PI / 180;
    double[3][3] rot1, rot2, rot;
    foreach(i; 0..3)
        foreach(j; 0..3)
            rot1[i][j] = rot2[i][j] = i==j ? 1.0 : 0.0;
    rot1[0][0] = cos(alpha); // xz rotation around y
    rot1[0][2] = -sin(alpha);
    rot1[2][0] = sin(alpha);
    rot1[2][2] = cos(alpha);
    rot2[1][1] = cos(beta); //yz rotation around x
    rot2[1][2] = -sin(beta);
    rot2[2][1] = sin(beta);
    rot2[2][2] = cos(beta);

    foreach(i;0..3)
        foreach(j;0..3)
            rot[i][j] = iota(3).map!(k => rot1[i][k] * rot2[k][j]).sum;
    return rot;
}

auto getOr(T, string mbr, V)(V def) {
    static if (__traits(hasMember, T, mbr)) return __traits(getMember, T, mbr);
    else return def;
}

class Renderer {
    abstract void drawSurface(ImageZ img, Params ps, bool big);
    abstract void drawPoints(ImageZ img, UV[] points, Params ps);
    abstract UV[] drawFloorRay(ImageZ img, Params ps, const double u0, const double v0);
    abstract UV walk(UV pos, ref double heading, double dt, double dist);
    abstract string name();
    abstract void setDistance();
    abstract void ensureCorrectPos(Params ps);
}

class Render(World) : Renderer {
    enum poleSingularity = __traits(hasMember, World, "singularity");
    alias Surf = Surface!(World.equation, poleSingularity);

    override string name() { return World.name; }
    override void setDistance() { globalDistance = World.distance; }
    override void ensureCorrectPos(Params ps) {
        static if (__traits(hasMember, World, "ensureCorrectPos"))
            World.ensureCorrectPos(ps);
    }
    override void drawSurface(ImageZ img, Params ps, bool big) {
        int NU = getOr!(World, "NU")(2000);
        int NV = getOr!(World, "NV")(1000); // number of mesh points
        if (big) {
            NV *= 2; NU = NU*3/2;
        }
        const W2 = img.W / 2, H2 = img.H / 2;
        const double R = globalR;
        const double zCenter = R * globalDistance, zScreen = img.W / 1.2;
        img.clear();
        double[3][3] rot = rotMatrix(ps);

        foreach(iv; 1..NV-1) {
            static if (__traits(hasMember, World, "pickV"))
            double v = World.pickV(iv, NV);
            else
            double v = (iv - NV/2) * PI / NV;

            foreach(iu; 0.. NU) {
                double u = cast(double)iu * 2.0*PI / NU;
                double x0,y0,z0;
                Surf.embedIn3D(u,v, x0,y0,z0);
                double x = rot[0][0]*x0 + rot[0][1]*y0 + rot[0][2]*z0;
                double y = rot[1][0]*x0 + rot[1][1]*y0 + rot[1][2]*z0;
                double z = rot[2][0]*x0 + rot[2][1]*y0 + rot[2][2]*z0;

                // project on 'screen' plane
                double z1 = zCenter + z;
                if (z1 < zScreen) continue;
                double k = zScreen / z1;
                double px = x * k, py = y * k;

                uint clr = World.color(u, v);
                int sx = W2 + cast(int)(px);
                int sy = H2 - cast(int)(py);
                img.putPixelZ(sx, sy, clr, z1);
            }
        }
    }

    ImageZ[] wrkImages;

    void ensureImagesCreated(size_t num, int w, int h) {
        wrkImages.length = num;
        foreach(ref im; wrkImages)
            if (im is null) im = new ImageZ(w,h);
    }

    override void drawPoints(ImageZ img, UV[] points, Params ps) {
        const W2 = img.W / 2, H2 = img.H / 2;
        const double zCenter = globalR * globalDistance, zScreen = img.W / 1.2;
        double[3][3] rot = rotMatrix(ps);

        bool toScr(UV p, ref int sx, ref int sy, ref double z) {
            double x0,y0,z0;
            Surf.embedIn3D(p.u, p.v, x0,y0,z0);
            double x = rot[0][0]*x0 + rot[0][1]*y0 + rot[0][2]*z0;
            double y = rot[1][0]*x0 + rot[1][1]*y0 + rot[1][2]*z0;
                   z = rot[2][0]*x0 + rot[2][1]*y0 + rot[2][2]*z0;
            double z1 = zCenter + z;
            if (z1 < zScreen) return false;
            double k = zScreen / z1;
            double px = x * k, py = y * k;
            sx = W2 + cast(int)(px);
            sy = H2 - cast(int)(py);
            return true;
        }
        int sx, sy;
        double z;
        if (toScr(ps.pos, sx, sy, z))
            foreach(dy; -1..2) // draw our position
                foreach(dx; -1..2)
                    img.putPixel(sx+dx, sy+dy, 0xFFFFFF);

        foreach(p; points) {
            if (!toScr(p, sx, sy, z)) continue;
            uint clr =  (128 + cast(int)(World.pointColor(z) * 120)) << 16;
            img.putPixel(sx, sy, clr);
        }
    }

    int[] whoseCol;

    override UV[] drawFloorRay(ImageZ img, Params ps, const double u0, const double v0) {
        import std.range : iota;
        import std.parallelism;
        const int W = img.W, H = img.H;
        const double screenZ = W / 1.2, camY = H*5/8;
        const int SkyH = H*3/8, FloorH = H*5/8;
        auto paths = appender!(UV[]);
        const bool dyndt = ps.dyndt, walls = ps.walls;
        const int minsy = cast(int)(FloorH / ps.range);
        ensureImagesCreated(taskPool.size + 1, W, H);

        if (walls) {
            foreach(im; wrkImages[].parallel)
                im.fillRows(0, SkyH + minsy, 0x505090);
        } else
            img.fillRows(0, SkyH + minsy, 0x505090);

        auto pathApps = iota(taskPool.size + 1).amap!(i => appender!(UV[]));
        whoseCol.length = W;
        const ray_sx = (-W/2) & 63;

        foreach(sx; iota(-W/2, W/2).parallel) {
            auto idx = taskPool.workerIndex;
            whoseCol[sx+W/2] = cast(int)idx;
            double dt = ps.dt;
            double angle = atan(sx / screenZ) + ps.heading*PI/180;
            const double scrDistAlongRay = sqrt(sx*sx + screenZ*screenZ);
            double s = 0, du,dv;
            Surf.courseVector(u0, v0, angle, du, dv);

            Dbl4 state;
            state.data[0] = u0;
            state.data[1] = du;
            state.data[2] = v0;
            state.data[3] = dv;

            int sy = FloorH, iters = 0;
            double nextDist = scrDistAlongRay;
            while(sy > minsy && iters < 5000) {
                const double u = state.data[0], v = state.data[2];

                if (walls) {
                    auto clr = World.color(u, v);
                    if (clr == wallColor) {
                        int totalh = cast(int) (H * scrDistAlongRay / s);
                        int h = min(totalh, H);
                        int y0 = SkyH - h*3/8;
                        assert(y0 >= 0);
                        int tx = cast(int)((u*cos(v) + v)*500) & 127;
                        double ky = 128.0 / totalh;
                        double ty0 = h==totalh ? 0 :  (totalh*3/8 - SkyH) * ky;
                        foreach(i; 0..h) {
                            int ty = cast(int) (ky*i + ty0);
                            assert(ty >= 0);
                            assert(ty < 128);
                            clr = wallTexture[tx*128 + ty];
                            wrkImages[idx].putPixel(sx + W/2, y0 + i, clr);
                        }
                        sy = minsy;
                        break;
                    }
                }

                if (s >= nextDist) {
                    auto clr = World.color(u, v);
                    wrkImages[idx].putPixel(sx + W/2, sy + SkyH-1, clr);
                    sy--;
                    nextDist = scrDistAlongRay * camY / sy;
                    if ((sx & 63)==ray_sx) pathApps[idx] ~= UV(u,v);
                    if (dyndt) dt = sqrt(nextDist - s);
                }
                double ds = Surf.vlen(u, v, state.data[1] * dt, state.data[3] * dt);
                state = evolveRK!(Surf.geodesicStep)(state, dt);
                s += ds;
                iters++;
            }
            while(sy > minsy) { // not full vertical line is drawn yet
                static if (poleSingularity)
                    wrkImages[idx].repeatPixel(sx + W/2, sy + SkyH-1);
                else
                    wrkImages[idx].putPixel(sx + W/2, sy + SkyH-1, 0);
                sy--;
            }
        }// for sx

        int y0 = walls ? 0 : SkyH+minsy;
        foreach(sx; 0..W) {
            auto idx = whoseCol[sx];
            foreach(y; y0 .. H) {
                auto j = y*W+sx;
                img.colors[j] = wrkImages[idx].colors[j];
            }
        }

        return pathApps.map!(a => a.data).join;
    }

    override UV walk(UV pos, ref double heading, double dt, double dist) {
        double angle = heading*PI/180;
        double s = 0, du,dv;
        Surf.courseVector(pos.u, pos.v, angle, du, dv);
        Dbl4 state;
        state.data[0] = pos.u;
        state.data[1] = du;
        state.data[2] = pos.v;
        state.data[3] = dv;
        int iters = 0;
        while(s < dist && iters < 1000) {
            double u = state.data[0], v = state.data[2];
            du = state.data[1] * dt; dv = state.data[3] * dt;
            double ds = Surf.vlen(u, v, du, dv);
            state = evolveRK!(Surf.geodesicStep)(state, dt);
            s += ds;
            iters++;
        }
        double u = state.data[0], v = state.data[2];
        double a = Surf.courseAngle(u,v, state.data[1], state.data[3]);
        heading = a*180/PI;
        while(heading < 0) heading += 360;
        while(heading >= 360) heading -= 360;
        return UV(u, v);
    }
}//Render
