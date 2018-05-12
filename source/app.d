import std.stdio, dlangui, dlangui.widgets.metadata, std.math, std.conv, std.array, std.format;
import symbolic;
mixin APP_ENTRY_POINT;

__gshared uint[] wallTexture;

class DrawingBoard : ImageWidget {
	ColorDrawBuf cdbuf;
    int W = 850, H = 600; //current size

	this() {
		cdbuf = new ColorDrawBuf(W, H);
		Ref!DrawBuf r = cdbuf;
		this.drawable = new ImageDrawable(r);
		this.margins = 0;// Rect(0, 10, 0, 10);//10
        fillBytes(0xc0c0c0);
        initWall();
	}

	override void measure(int parentWidth, int parentHeight) {
        measuredContent(parentWidth, parentHeight, W, H);
    }

    override void onDraw(DrawBuf buf) {
        if (visibility != Visibility.Visible)
            return;
        super.onDraw(buf);
    }

	void fillBytes(uint q) {
		foreach(y; 0..cdbuf.height) {
			uint* ptr = cdbuf.scanLine(y);
			foreach(x; 0..cdbuf.width)
				ptr[x] = q;// ((x + q) ^ y) & 255;
		}
	}

    void drawImgZ(ImageZ img) {
        auto w = img.W;
        foreach(y; 0..img.H) {
            uint* ptr = cdbuf.scanLine(y);
            ptr[0..w] = img.colors[y*w .. y*w+w];
        }
        invalidate();
    }

    void drawImgAt(ImageZ img, int dx, int dy, int sx, int sy, int szx, int szy) {
        assert(dx >= 0 && dy >= 0 && dx + szx <= W && dy + szy <= H);
        auto w = img.W;
        foreach(y; 0..szy) {
            uint* ptr = cdbuf.scanLine(dy + y);
            ptr[dx .. dx + szx] = img.colors[(sy+y)*w + sx .. (sy+y)*w + sx + szx];
        }
    }

    void initWall() {
        import std.random : uniform;
        wallTexture.length = 128*128;
        foreach(ref x; wallTexture) {
            int r = uniform(0, 20), g = uniform(0, 20), b = uniform(0, 20);
            x = 0xB0503F + (r<<16) + (g<<8) + b;
        }
        foreach(row; 0..8) {
            foreach(y; 0..2)
            foreach(x; 0..128) {
                int r = uniform(0, 20), g = uniform(0, 20), b = uniform(0, 20);
                wallTexture[row*16+y+x*128] = 0xCDC5BA  + (r<<16) + (g<<8) + b;
            }
            int dx = (row&1)*16;
            foreach(n; 0..4)
            foreach(y; 0..16)
                foreach(x; 0..2) {
                    int r = uniform(0, 20), g = uniform(0, 20), b = uniform(0, 20);
                    wallTexture[row*16+y+ (dx + n*32 + x)*128] = 0xCDC5BA  + (r<<16) + (g<<8) + b;
                }
        }
        foreach(y; 0..128) {
            uint* ptr = cdbuf.scanLine(y);
            ptr[0..128] = wallTexture[y*128 .. y*128+128];
        }
    }
}

mixin(registerWidgets!("__gshared static this", DrawingBoard));

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
__gshared double globalR = 1000.0;

string[] genCode(alias surfaceEquation)() {
    Expr[] X = surfaceEquation();
    string embed = format("x = %s; y = %s; z = %s;", X[0].code, X[1].code, X[2].code);

    auto Xu = X.diff("u"); // basis vectors
    auto Xv = X.diff("v");
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
                    tmp[l] = div( mul(ginv[k][l], new Add([J, I, L]).simp), two);
                }
                C[k][i][j] = add(tmp[0], tmp[1]);
            }
    // res.data[1] = -C[0][0][0]*du*du - 2 * C[0][0][1]*du*dv - C[0][1][1]*dv*dv; //u''
    // res.data[3] = -C[1][0][0]*du*du - 2 * C[1][0][1]*du*dv - C[1][1][1]*dv*dv; //v''
    string u_accel = format("res.data[1] = -%s*du*du - 2 * %s*du*dv - %s*dv*dv;\n", 
                            C[0][0][0].code, C[0][0][1].code, C[0][1][1].code).txtSimp; 
    string v_accel = format("res.data[3] = -%s*du*du - 2 * %s*du*dv - %s*dv*dv;\n", 
                            C[1][0][0].code, C[1][0][1].code, C[1][1][1].code).txtSimp; 

    Expr du = new Var("du"), dv = new Var("dv");
    Expr mag = add(mul(g[0][0], du, du), mul(two, g[0][1], du, dv), mul(g[1][1], dv,dv));
    string vlen = format("return sqrt(%s);", mag.code).txtSimp;                            
    return [embed, u_accel, v_accel, vlen, g[0][0].code, g[0][1].code, g[1][1].code];
    //      0      1        2         3      4            5              6
}

class Surface(alias surfaceEquation) {
    enum codes = genCode!surfaceEquation();
    
    static void embedIn3D(double u, double v, ref double x, ref double y, ref double z) {
        pragma(msg, codes[0]);
        const double R = globalR;
        const double cos_v = cos(v);
        const double sin_v = sin(v);
        const double cos_u = cos(u);
        const double sin_u = sin(u);
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
        if (abs(cos_v) < 0.0000001) { // handle poles, poorly
            res.data[1] = du;
            res.data[3] = dv;
            return res;
        }
        const double sin_v = sin(v);
        const double sin_u = sin(u), cos_u = cos(u);
        pragma(msg, codes[1]);
        pragma(msg, codes[2]);
        mixin(codes[1]);
        mixin(codes[2]);
        return res;
    }

    static uint color(double u, double v) {
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

        // double l = (sin(v) + 1.0)*0.5;
        int c = ((iu ^ iv) & 31)*4 + 64;
        uint clr = (c << 16) + (c<<8) + c;
        int mu = iu % 20, mv = (iv + 180) % 20;
        if (mu >= 10 && mv >= 10 && mu < 16 && mv < 16)
            clr = wallColor; //squares
        else {
            if (mv >= 13 && mv < 16 && mu >= 3 && mu < 6)
                clr = (iu/3 * 2)<<16;
            else
            if (mv >= 3 && mv < 6 && mu >= 13 && mu < 16)
                clr = ((iv + 90)/3 * 4)<<8;
        }
        return clr;
    }

    static double vlen(double u, double v, double du, double dv) {
        pragma(msg, codes[3]);
        const double R = globalR;
        const double cos_v = cos(v), sin_v = sin(v);
        const double cos_u = cos(u), sin_u = sin(u);
        mixin(codes[3]);
    }

    static void courseVector(double u0, double v0, double angle, ref double du, ref double dv) {
        const double R = globalR;
        double cos_v = cos(v0), sin_v = sin(v0);
        double cos_u = cos(u0), sin_u = sin(u0);
        if (abs(cos_v) < 0.0000001) cos_v = 1.0; // nonsensical but at least won't crash
        enum u_code = format("du = sin(angle) / sqrt(%s);", codes[4]).txtSimp;
        enum v_code = format("dv = cos(angle) / sqrt(%s);", codes[6]).txtSimp;
        pragma(msg, u_code);
        pragma(msg, v_code);
        mixin(u_code);
        mixin(v_code);
        //du = sin(angle)/(cos_v*R); dv = cos(angle)/R;
    }

    static double courseAngle(double u0, double v0, double du, double dv) { // in radians
        // up (0, 1/sqrt(g22=R^2)) = (0, 1/R)
        // up*V = cos(a) * |V| = g12 * up_v * du + g22 * up_v * dv = R * dv
        // cos(a) = R*dv / |V|
        const double R = globalR;
        const double cos_v = cos(v0), sin_v = sin(v0);
        const double sin_u = sin(u0), cos_u = cos(u0);
        enum up_v_code = format("double up_v = 1/sqrt(%s);", codes[6]).txtSimp;
        enum product_code = format("double product = %s * up_v * du + %s * up_v * dv;", codes[5], codes[6]).txtSimp;
        pragma(msg, up_v_code);
        pragma(msg, product_code);
        mixin(up_v_code);
        mixin(product_code);
        double a = acos(product / vlen(u0, v0, du, dv));
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

Expr[] thingyEq() {
    auto R = new Var("R"), two = new Const("2");
    /*
    X = 2+cos(v)
    z = sin(v)
    x = X*cos(u) = (2+cos(v))*cos(u)
    y = X*sin(u) = (2+cos(v))*sin(u)
    */

    return [mul(R, add(two, new Cos("v")), new Cos("u")),
            mul(R, add(two, new Cos("v")), new Sin("u")),
            mul(R, new Sin("v")) ];
}

class Params {
    double heading, dt, range;
    UV pos;
    bool dyndt, walls;
}

class Renderer {
    abstract void drawSurface(ImageZ img);
    abstract void drawPoints(ImageZ img, UV[] points);
    abstract UV[] drawFloorRay(ImageZ img, Params ps, double u0, double v0);
    abstract UV walk(UV pos, ref double heading, double dt, double dist);
}

class Render(alias surfaceEq) : Renderer {
    alias Surf = Surface!surfaceEq;

    override void drawSurface(ImageZ img) {
        enum NU = 1500, NV = 700; // number of mesh points
        const W2 = img.W / 2, H2 = img.H / 2;
        const double R = globalR;
        const double zCenter = 9*R, zScreen = img.W / 1.2;
        foreach(iv; 1..NV-1) {
            double v = cast(double)(iv - NV/2) * 2* PI / NV;
            foreach(iu; 0.. NU) {
                double u = cast(double)iu * 2.0*PI / NU;
                double x,y,z;
                Surf.embedIn3D(u,v, x,y,z);

                // project on 'screen' plane
                double z1 = zCenter + z;
                double k = zScreen / z1;
                double px = x * k, py = y * k;

                uint clr = Surf.color(u, v);
                int sx = W2 + cast(int)(px);
                int sy = H2 - cast(int)(py);
                img.putPixelZ(sx, sy, clr, z1);
            }
        }
    }

    override void drawPoints(ImageZ img, UV[] points) {
        const W2 = img.W / 2, H2 = img.H / 2;
        const double zCenter = 9 * globalR, zScreen = img.W / 1.2;
        foreach(p; points) {
            double x,y,z;
            Surf.embedIn3D(p.u, p.v, x,y,z);
            double z1 = zCenter + z;
            double k = zScreen / z1;
            double px = x * k, py = y * k;
            uint clr = (128 - cast(int)(z/globalR * 100)) << 16;
            int sx = W2 + cast(int)(px);
            int sy = H2 - cast(int)(py);
            img.putPixel(sx, sy, clr);
        }
    }


    override UV[] drawFloorRay(ImageZ img, Params ps, double u0, double v0) {
        import std.range : iota;
        const int W = img.W, H = img.H;
        double screenZ = W / 1.2, camY = H/2;
        const double x0 = u0, z0 = v0;
        auto paths = appender!(UV[]);
        const bool dyndt = ps.dyndt, walls = ps.walls;
        const int minsy = cast(int)(H/2 / ps.range);
        img.fillRows(0, H/2+minsy, 0x505090);

        foreach(sx; iota(-W/2, W/2)) {
            double dt = ps.dt;
            double angle = atan(sx / screenZ) + ps.heading*PI/180;
            double scrDistAlongRay = sqrt(sx*sx + screenZ*screenZ);
            double s = 0, dx,dz;
            Surf.courseVector(x0, z0, angle, dx, dz);

            Dbl4 state;
            state.data[0] = x0;
            state.data[1] = dx;
            state.data[2] = z0;
            state.data[3] = dz;

            int sy = H/2, iters = 0;
            double nextDist = scrDistAlongRay;
            while(sy > minsy && iters < 5000) {
                double x = state.data[0], z = state.data[2];

                if (walls) {
                    auto clr = Surf.color(x, z);
                    if (clr == wallColor) {
                        int totalh = cast(int) (H * scrDistAlongRay / s);
                        int h = totalh <= H-2 ? totalh : H-2;
                        int y0 = H/2-1-h/2;
                        int tx = cast(int)((x*cos(z) + z)*500) & 127;
                        double ky = 128.0 / totalh;
                        double ty0 = h==totalh ? 0 :  (totalh-h)/2 * ky;
                        foreach(i; 0..h) {
                            int ty = cast(int) (ky*i + ty0);
                            clr = wallTexture[tx*128 + ty];
                            img.putPixel(sx + W/2, y0 + i, clr);
                        }
                        break;
                    }
                }

                if (s >= nextDist) {
                    auto clr = Surf.color(x, z);
                    img.putPixel(sx + W/2, sy + H/2-1, clr);
                    sy--;
                    nextDist = scrDistAlongRay * camY / sy;
                    if ((sx & 63)==0) paths ~= UV(x,z);
                    if (dyndt) dt = sqrt(nextDist - s);
                }
                double du = state.data[1] * dt, dv = state.data[3] * dt;
                double ds = Surf.vlen(x, z, du, dv);
                state = evolveRK!(Surf.geodesicStep)(state, dt);
                s += ds;
                iters++;
            }
        }// for sx
        return paths.data;
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
}

extern (C) int UIAppMain(string[] args) {
    import dlangui.core.logger;
    int w= 880, h = 700;
   	Log.setLogLevel( dlangui.core.logger.LogLevel.Error );

    version(Windows) {
        w = w.pixelsToPoints; h = h.pixelsToPoints;
    }
    Window window = Platform.instance.createWindow("Curved", null, 1, w, h);
    window.mainWidget = parseML(q{
        VerticalLayout {
            backgroundColor: 0xc0c0c0
            // margins: 10pt
            padding: 10pt
            layoutWidth: 750
            HorizontalLayout {
                Button {text: "Render"; id: "btnRender"}
                TextWidget {text: "Heading:" }
                EditLine { text: "90"; id: "heading"; layoutWidth: 50}
                TextWidget {text: "dt:" }
                EditLine { text: "1"; id: "dt"; layoutWidth: 50}
                CheckBox { text: "Dyn.dt"; id: "dyndt"}
                CheckBox { text: "Walls"; id: "walls"}
                TextWidget {text: "Range:" }
                EditLine { text: "6"; id: "range"; layoutWidth: 70}
                TextWidget {text: "R:" }
                EditLine { text: "1000"; id: "R"; layoutWidth: 70}                
                TextWidget {text:""; id:"out"}
            }
            DrawingBoard {id: "pic"}
        }
    });

    auto edHeading = window.mainWidget.childById!EditLine("heading");
    auto edDt = window.mainWidget.childById!EditLine("dt");
    auto txtOut = window.mainWidget.childById!TextWidget("out");
    auto cbDDT = window.mainWidget.childById!CheckBox("dyndt");
    auto cbWalls = window.mainWidget.childById!CheckBox("walls");
    auto edRange = window.mainWidget.childById!EditLine("range");
    auto edR = window.mainWidget.childById!EditLine("R");
    auto imgSphere = new ImageZ(512,512);
    auto img = new ImageZ(512,512);
    auto imgFloor = new ImageZ(512,512);
    auto pic = window.mainWidget.childById!DrawingBoard("pic");
    auto ps = new Params();
    ps.pos.u = 4.7; ps.pos.v = 0.2;
    ps.range = 6.0;
    Renderer rend = new Render!thingyEq;
    rend.drawSurface(imgSphere);

    void render() {
        import std.datetime.stopwatch : StopWatch, AutoStart;
        auto sw = StopWatch(AutoStart.yes);
        ps.heading = edHeading.text.to!double;
        ps.dt = edDt.text.to!double;
        ps.dyndt = cbDDT.checked;
        ps.walls = cbWalls.checked;
        auto rng = edRange.text.to!double;
        if (rng >= 1 && rng < 10) ps.range = rng;
        auto r = edR.text.to!double;
        if (r >= 100 && r <= 1000000) globalR = r; 
        img.copyFrom(imgSphere);
        auto points = rend.drawFloorRay(imgFloor, ps, ps.pos.u, ps.pos.v);
        rend.drawPoints(img, points);
        pic.drawImgAt(img, 0,0, 100,100, 320,320);
        pic.drawImgAt(imgFloor, 320,0, 0,0, 512,512);
        txtOut.text = format("t=%s"d, sw.peek.total!"msecs");
    }

    window.mainWidget.childById!Button("btnRender").click = delegate(Widget w) {
        render();
        return true;
    };

    window.mainWidget.keyEvent = delegate(Widget wt, KeyEvent e) {
        //writeln("main keyEvent ", e);
        if (e.action==KeyAction.KeyDown) {
            switch(e.keyCode) {
                case KeyCode.KEY_A: ps.heading -= 10; edHeading.text = ps.heading.to!dstring; break;
                case KeyCode.KEY_D: ps.heading += 10; edHeading.text = ps.heading.to!dstring; break;
                case KeyCode.KEY_W: 
                    auto hdn = ps.heading;
                    ps.pos = rend.walk(ps.pos, hdn, ps.dt, 50); 
                    ps.heading = hdn;
                    edHeading.text = format("%.3g"d, hdn);
                    break;
                case KeyCode.KEY_S: 
                    auto hdn = ps.heading + 180;
                    ps.pos = rend.walk(ps.pos, hdn, ps.dt, 50); 
                    ps.heading = hdn - 180;
                    while(ps.heading < 0) ps.heading += 360;
                    while(ps.heading >= 360) ps.heading -= 360;
                    edHeading.text = format("%.3g"d, ps.heading);
                    break;
                default: return false;
            }
            render();
            window.invalidate();
            return true;
        }
        return false;
    };

    window.show();
    return Platform.instance.enterMessageLoop();
}
