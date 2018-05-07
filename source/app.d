import std.stdio, dlangui, dlangui.widgets.metadata, std.math, std.conv, std.array, std.format;

mixin APP_ENTRY_POINT;

class MotionPicture : ImageWidget {
	ColorDrawBuf cdbuf;
    int W = 700, H = 600; //current size

	this() {
		cdbuf = new ColorDrawBuf(W, H);
		Ref!DrawBuf r = cdbuf;
		this.drawable = new ImageDrawable(r);
		this.margins = 0;// Rect(0, 10, 0, 10);//10
        fillBytes(0x404040);
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
}

mixin(registerWidgets!("__gshared static this", MotionPicture));

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
}

void drawSphere(ImageZ img) {
    // R = 1   u - longitude,  v - latitude
    // x = cos(u)*cos(v)
    // z = sin(u)*cos(v)   // swapped y and z
    // y = sin(v)
    enum NU = 1000, NV = 500; // number of mesh points
    const W2 = img.W / 2, H2 = img.H / 2;
    const double R = Sphere.R;
    const double zCenter = 3*R, zScreen = img.W / 1.2;
    foreach(iv; 1..NV-1) {
        double v = cast(double)(iv - NV/2) * PI / NV;
        foreach(iu; 0.. NU) {
            double u = cast(double)iu * 2.0*PI / NU;
            double x,y,z;
            Sphere.embedIn3D(u,v, x,y,z);

            // project on 'screen' plane
            double z1 = zCenter + z;
            double k = zScreen / z1;
            double px = x * k, py = y * k;

            //int ux = cast(int)(u * 100), vx = cast(int)(v * 100) + 1024;
            uint clr = Sphere.color(u, v);
            int sx = W2 + cast(int)(px);
            int sy = H2 - cast(int)(py);
            img.putPixelZ(sx, sy, clr, z1);
            // f.writefln("iu=%s iv=%s x,y,z=%s,%s,%s px,py=%s,%s sx,sy=%s,%s, clr=%s",
            //             iu, iv,  x,y,z,  px,py,   sx,sy, clr & 255);

            // int c = /*((iu ^ iv) & 127) +*/ ((cast(int) (v/PI*100) + 64)&127);
            // f.writefln("u=%s v=%s v/PI*100=%s c=%s", u,v, v/PI*100, c);
        }
    }
}

UV[] makeOneSphereGeodesic() {
    Dbl4 state;
    state.data = [-1.5, 0.1, 0.0, 0.1]; // one step ds = dt * speed
    double t = 0, dt = 2*PI / 360.0, s = 0.0;
    auto points = appender!(UV[]);

    foreach(i; 0..3600) {
        auto u = state.data[0], v = state.data[2];
        double du = state.data[1] * dt, dv = state.data[3] * dt;
        double ds = Sphere.vlen(u,v,du,dv);
        points ~= UV(u,v);
        state = evolveRK!(Sphere.geodesicStep)(state, dt);
        t += dt;
        s += ds;
        if (s > 6280.0) break;
    }
    return points.data;
}

struct UV { double u,v; }

void drawPoints(S)(ImageZ img, UV[] points) {
    const W2 = img.W / 2, H2 = img.H / 2;
    const double zCenter = 3 * Sphere.R, zScreen = img.W / 1.2;
    foreach(p; points) {
        double x,y,z;
        S.embedIn3D(p.u, p.v, x,y,z);
        double z1 = zCenter + z;
        double k = zScreen / z1;
        double px = x * k, py = y * k;
        uint clr = (128 - cast(int)(z/Sphere.R * 100)) << 16;
        int sx = W2 + cast(int)(px);
        int sy = H2 - cast(int)(py);
        img.putPixel(sx, sy, clr);
    }
}

// time-dependent:  f'(t,y) = drv(t,y)
T evolveRKt(alias drv, T)(double t, T state, double dt) {
    auto a = drv(t, state);
    auto b = drv(t + dt/2,  state + a * (dt/2));
    auto c = drv(t + dt/2,  state + b * (dt/2));
    auto d = drv(t + dt, state + c*dt);
    return state + a*(dt/6) + b*(dt/3) + c*(dt/3) + d*(dt/6);
}

// f'(t,y) = drv(y)
T evolveRK(alias drv, T)(T state, double dt) {
    auto a = drv(state);
    auto b = drv(state + a * (dt/2));
    auto c = drv(state + b * (dt/2));
    auto d = drv(state + c*dt);
    return state + a*(dt/6) + b*(dt/3) + c*(dt/3) + d*(dt/6);
}

double dsqrt(double x, double y) {
    return 1.0/(2*sqrt(x));
}

void testRK() {
    double s = 1.0;
    double dt = 1.0/128, t = 1.0;
    foreach(n; 0..3*128+1) {
        if (n % 128 == 0) {
            writefln("x=%s y=%s", t, s);
        }
        s = evolveRKt!dsqrt(t, s, dt);
        t += dt;
    }
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

Dbl4 geod(Dbl4 y) {
    //y : u, u', v, v'
    Dbl4 res;
    auto du = res.data[0] = y.data[1]; // u' = u'
    auto dv = res.data[2] = y.data[3]; // v' = v'
    static double[2][2][2] C;
    //   [[(0, -tan(v)), (0,0)],[(sin(v)*cos(v),0),(0, 0)]])
    //    C = C.subs({u:y0,v:y2})
    double v = y.data[2];
    C[0][0][0] = 0; // sphere
    C[0][0][1] = -tan(v);
    C[0][1][0] = 0;
    C[0][1][1] = 0;
    C[1][0][0] = sin(v)*cos(v);
    C[1][0][1] = 0;
    C[1][1][0] = 0;
    C[1][1][1] = 0;

    res.data[1] = -C[0][0][0]*du*du - 2 * C[0][0][1]*du*dv - C[0][1][1]*dv*dv; //u''
    res.data[3] = -C[1][0][0]*du*du - 2 * C[1][0][1]*du*dv - C[1][1][1]*dv*dv; //v''

    return res;
}

class FlatPlane {
    static Dbl4 geodesicStep(Dbl4 y) {
        //y : u, u', v, v'
        Dbl4 res;
        auto du = res.data[0] = y.data[1]; // u' = u'
        auto dv = res.data[2] = y.data[3]; // v' = v'
        res.data[1] = 0; // u'' = 0
        res.data[3] = 0;
        return res;
    }

    static uint color(double u, double v) {
        int iu = cast(int)u, iv = cast(int)v;
        int c = ((iu ^ iv) & 255) | 128;
        uint clr = (c << 16) + (c << 8) + c;
        return clr;
    }

    static double vlen(double u, double v, double du, double dv) {
        return sqrt(du*du + dv*dv);
    }

    static void courseVector(double u0, double v0, double angle, ref double du, ref double dv) {
        du = sin(angle); dv = cos(angle);
    }
}//FlatPlane

class Sphere {
/*
g = (cos(v)^^2  0
     0          1.0)
    C[0][0][1] = -tan(v);
    C[1][0][0] = sin(v)*cos(v);
    C others = 0
    res.data[1] = -C[0][0][0]*du*du - 2 * C[0][0][1]*du*dv - C[0][1][1]*dv*dv; //u''
    res.data[3] = -C[1][0][0]*du*du - 2 * C[1][0][1]*du*dv - C[1][1][1]*dv*dv; //v''

*/
    static double R = 1000.0;
    static void embedIn3D(double u, double v, ref double x, ref double y, ref double z) {
        x = R*cos(u)*cos(v);
        z = R*sin(u)*cos(v);
        y = R*sin(v);// - R*0.8;
    }

    static Dbl4 geodesicStep(Dbl4 y) {
        //y : u, u', v, v'
        Dbl4 res;
        auto du = res.data[0] = y.data[1]; // u' = u'
        auto dv = res.data[2] = y.data[3]; // v' = v'
        double v = y.data[2];
        double cosv = cos(v);
        if (abs(cosv) < 0.0000001) {
            res.data[1] = du;
            res.data[3] = dv;
            return res;
        }
        double sinv = sin(v);
        double C001 = -sinv/cosv, C100 = sinv*cosv;
        res.data[1] = -2 * C001*du*dv; //u''
        res.data[3] = -C100*du*du; //v''
        return res;
    }

    static uint color(double u, double v) {
        int iu = cast(int)(u/PI*180), iv = cast(int)(v/PI*180);
        if (iu >= 360) {
            do iu -= 360; while (iu >= 360);
        } else if (iu < 0) {
            do iu += 360; while (iu < 0);
        }

        if (iv >= 90) {
            iv = 180 - iv;
        } else if (iv < -90) {
            iv = -180 - iv;
        }

        double l = (sin(v) + 1.0)*0.5;
        int c = cast(int)(((iu ^ iv) & 255)*l);
        uint clr = /*(c << 16) + (c << 8) +*/ c;
        if (((iu/10)&1)==1 && (((iv+1000)/10)&1)==1) clr = 0xF08000; //red squares
        return clr;
    }

    static double vlen(double u, double v, double du, double dv) {
        double g11 = cos(v)^^2;
        return sqrt(g11*du*du + dv*dv)*R;
    }

    static void courseVector(double u0, double v0, double angle, ref double du, ref double dv) {
        double cosv0 = cos(v0);
        if (abs(cosv0) < 0.0000001) cosv0 = 1.0; // nonsensical but at least won't crash
        du = sin(angle)/(cosv0*R); dv = cos(angle)/R;
    }

    static double courseAngle(double u0, double v0, double du, double dv) { // in radians
        // up (0, 1/sqrt(g22=R^2)) = (0, 1/R)
        // up*V = cos(a) * |V| = g12 * up_v * du + g22 * up_v * dv = R * dv
        // cos(a) = R*dv / |V|
        double a = acos(R*dv / vlen(u0, v0, du, dv));
        if (du < 0) a = 2*PI - a;
        return a;
    }
}//Sphere

void drawFloorDirect(ImageZ img, Params ps) {
    const int W = img.W, H = img.H;
    double screenZ = W / 1.2, camY = H/2;
    foreach(sx; -W/2 .. W/2) {
        foreach(sy; 40.. H/2) {
            double k = camY / sy;
            double x = sx * k;
            double z = screenZ * k;
            int ix = cast(int)x, iz = cast(int)z;
            int c = ((ix ^ iz) & 255) | 128;
            uint clr = (c << 16) + (c << 8) + c;
            img.putPixel(sx + W/2, sy + H/2, clr);
        }
    }
}

UV[] drawFloorRay(S)(ImageZ img, Params ps, S space, double u0, double v0) {
    import std.range : iota;
    const int W = img.W, H = img.H;
    double screenZ = W / 1.2, camY = H/2;
    auto f = File("rays.txt", "wt");
    const double x0 = u0, z0 = v0;
    auto paths = appender!(UV[]);
    const int pixelWidth = 1;
    const bool dyndt = ps.dyndt;

    foreach(sx; iota(-W/2, W/2, pixelWidth)) {
        double dt = ps.dt;
        double angle = atan(sx / screenZ) + ps.heading*PI/180;
        double scrDistAlongRay = sqrt(sx*sx + screenZ*screenZ);
        double s = 0, t=0, dx,dz;
        space.courseVector(x0, z0, angle, dx, dz);

        Dbl4 state;
        state.data[0] = x0;
        state.data[1] = dx;
        state.data[2] = z0;
        state.data[3] = dz;

        const int minsy = 40;//H/2/8;
        int sy = H/2, iters = 0;
        double nextDist = scrDistAlongRay;
        bool chat = sx==0;
        int lastIter = 0;
        while(sy > minsy && iters < 5000) {
            double x = state.data[0], z = state.data[2];
            if (s >= nextDist) {
                auto clr = space.color(x, z);
                foreach(d; 0..pixelWidth)
                    img.putPixel(sx + W/2 + d, sy, clr);
                sy--;
                nextDist = scrDistAlongRay * camY / sy;
                if (chat) {
                    f.writefln("n=%s s=%s nestStep=%s", iters - lastIter, s, nextDist - s);
                    lastIter = iters;
                }
                if ((sx & 63)==0) paths ~= UV(x,z);
                if (dyndt) dt = sqrt(nextDist - s);
            }
            double du = state.data[1] * dt, dv = state.data[3] * dt;
            double ds = space.vlen(x, z, du, dv);
            state = evolveRK!(space.geodesicStep)(state, dt);
            t += dt;
            s += ds;
            iters++;
        }
        //f.writefln("sx=%s n=%s", sx, t / dt);
    }// for sx
    return paths.data;
}

UV walk(S)(UV pos, ref double heading, double dt, double dist) {
    double angle = heading*PI/180;
    double s = 0, du,dv;
    S.courseVector(pos.u, pos.v, angle, du, dv);
    Dbl4 state;
    state.data[0] = pos.u;
    state.data[1] = du;
    state.data[2] = pos.v;
    state.data[3] = dv;
    int iters = 0;
    while(s < dist && iters < 1000) {
        double u = state.data[0], v = state.data[2];
        du = state.data[1] * dt; dv = state.data[3] * dt;
        double ds = S.vlen(u, v, du, dv);
        state = evolveRK!(S.geodesicStep)(state, dt);
        s += ds;
        iters++;
    }
    double u = state.data[0], v = state.data[2];
    double a = S.courseAngle(u,v, state.data[1], state.data[3]);
    heading = a*180/PI;
    return UV(u, v);
}

class Params {
    double heading, dt;
    UV pos;
    bool dyndt;
}

extern (C) int UIAppMain(string[] args) {
    import dlangui.core.logger;
    int w= 800, h = 700;
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
                TextWidget {text:""; id:"out"}
            }
            MotionPicture {id: "pic"}
        }
    });

    // testRK();
    // Vec!(double, 4) v;
    // v.data[] = [1.0,2,3,4];
    // auto v2 = v * 1.1 + v;
    // writeln(v2.data);
    auto edHeading = window.mainWidget.childById!EditLine("heading");
    auto edDt = window.mainWidget.childById!EditLine("dt");
    auto txtOut = window.mainWidget.childById!TextWidget("out");
    auto cbDDT = window.mainWidget.childById!CheckBox("dyndt");

    auto imgSphere = new ImageZ(512,512);
    auto img = new ImageZ(512,512);
    auto imgFloor = new ImageZ(512,512);
    auto pic = window.mainWidget.childById!MotionPicture("pic");
    auto ps = new Params();
    ps.pos.u = 4.7; ps.pos.v = 0;
    auto plane = new FlatPlane();
    auto sphere = new Sphere();
    drawSphere(imgSphere);

    void render() {
        import std.datetime.stopwatch : StopWatch, AutoStart;
        //writeln("render");
        auto sw = StopWatch(AutoStart.yes);
        ps.heading = edHeading.text.to!double;
        ps.dt = edDt.text.to!double;
        ps.dyndt = cbDDT.checked;
        //drawSphere(imgSphere);
        img.copyFrom(imgSphere);
        //drawOneSphereGeodesic(img);
        //auto points = makeOneSphereGeodesic();
        //drawFloorDirect(imgZ, ps);
        auto points = drawFloorRay(imgFloor, ps, sphere, ps.pos.u, ps.pos.v);
        drawPoints!Sphere(img, points);
        pic.drawImgAt(img, 0,0, 100,100, 320,320);
        pic.drawImgAt(imgFloor, 0,320, 0,0, 512,256);
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
                    ps.pos = walk!Sphere(ps.pos, hdn, ps.dt, 50); 
                    ps.heading = hdn;
                    edHeading.text = format("%.3g"d, hdn);
                    break;
                case KeyCode.KEY_S: 
                    auto hdn = ps.heading + 180;
                    ps.pos = walk!Sphere(ps.pos, hdn, ps.dt, 50); 
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
