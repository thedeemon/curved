import std.stdio, dlangui, dlangui.widgets.metadata, std.math;

mixin APP_ENTRY_POINT;

class MotionPicture : ImageWidget {
	ColorDrawBuf cdbuf;
    int W = 700, H = 600; //current size

	this() {
		cdbuf = new ColorDrawBuf(W, H);
		Ref!DrawBuf r = cdbuf;
		this.drawable = new ImageDrawable(r);
		this.margins = 0;// Rect(0, 10, 0, 10);//10
        fillBytes(0);
	}

	override void measure(int parentWidth, int parentHeight) {
        /*DrawableRef img = drawable;
        int w = 0;
        int h = 0;
        if (!img.isNull) {
            w = showW;
            h = showH;
			if (fitImage) {
				auto mgs = margins;
				auto sz = imgSizeScaled(parentWidth, parentHeight - mgs.top - mgs.bottom);
				w = sz.x; h = sz.y + mgs.top + mgs.bottom;
			}
        }*/
        measuredContent(parentWidth, parentHeight, W, H);
    }

    override void onDraw(DrawBuf buf) {
        if (visibility != Visibility.Visible)
            return;
        super.onDraw(buf);
        // Rect rc = _pos;
        // applyMargins(rc);
		// auto saver = ClipRectSaver(buf, rc, alpha);
		// applyPadding(rc);
        // DrawableRef img = drawable;
        // if (!img.isNull) {
        //     //img.drawTo(buf, rc, state);
        // }
    }

	void fillBytes(int q) {
		foreach(y; 0..cdbuf.height) {
			uint* ptr = cdbuf.scanLine(y);
			foreach(x; 0..cdbuf.width)
				ptr[x] = ((x + q) ^ y) & 255;
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
    enum NU = 200, NV = 100;
    const W2 = img.W / 2, H2 = img.H / 2;
    const double zCenter = 2.0, zScreen = 1.0;

    foreach(iv; 1..NV-1) {
        double v = cast(double)(iv - NV/2) * PI / NV;
        foreach(iu; 0.. NU) {
            double u = cast(double)iu * 2.0*PI / NU;

            double x = cos(u)*cos(v);
            double z = sin(u)*cos(v);
            double y = sin(v);


            // project on 'screen' plane
            double z1 = zCenter + z;
            double k = zScreen / z1;
            double px = x * k, py = y * k;

            int ux = cast(int)(u * 100), vx = cast(int)(v * 100) + 1024;
            uint clr = ux ^ vx;
            int sx = W2 + cast(int)(px*W2);
            int sy = H2 - cast(int)(py*W2);
            img.putPixelZ(sx, sy, clr, z1);
        }
    }

    auto f = File("out.txt", "wt");
    Dbl4 state;
    state.data = [-1.5, 0.1, 0.0, 0.1]; // here speed = 0.141421;  one step ds = dt * speed
    double t = 0, dt = 2*PI / 360.0, s = 0.0;
    foreach(i; 0..3600) {
        auto u = state.data[0], v = state.data[2];

        double g11 = cos(v) ^^ 2;
        double g22 = 1.0;
        double du = state.data[1] * dt, dv = state.data[3] * dt;
        double ds = sqrt(g11*du*du + g22*dv*dv);

        double x = cos(u)*cos(v);
        double z = sin(u)*cos(v);
        double y = sin(v);


        // project on 'screen' plane
        double z1 = zCenter + z;
        double k = zScreen / z1;
        double px = x * k, py = y * k;

        int ux = cast(int)(u * 100), vx = cast(int)(v * 100) + 1024;
        uint clr = (128 - cast(int)(z*100)) << 16;
        int sx = W2 + cast(int)(px*W2);
        int sy = H2 - cast(int)(py*W2);
        img.putPixel(sx, sy, clr);

        if (i % 15 == 0) f.writefln("%s u=%s v=%s s=%s", i, u, v,s);
        state = evolveRK!geod(state, dt);
        t += dt;
        s += ds;
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

Dbl4 geodPlane(Dbl4 y) {
    //y : u, u', v, v'
    Dbl4 res;
    auto du = res.data[0] = y.data[1]; // u' = u'
    auto dv = res.data[2] = y.data[3]; // v' = v'
    res.data[1] = 0; // u'' = 0
    res.data[3] = 0;
    return res;    
}

void drawFloorDirect(ImageZ img) {
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

void drawFloorRay(ImageZ img) {
    const int W = img.W, H = img.H;
    double screenZ = W / 1.2, camY = H/2;
    const double dt = 0.5; 
    foreach(sx; -W/2 .. W/2) {
        double x = sx, z = screenZ;
        int lastsy = 999999; //target sy
        double k = camY / (camY-1);
        double dx = sx * k - x;
        double dz = screenZ * k - z;
        Dbl4 state;
        state.data[0] = x;
        state.data[1] = dx;
        state.data[2] = z;
        state.data[3] = dz;
        double s = sqrt(x*x + z*z);
        double scrDistAlongRay = s, t = 0;
        
        // |v0| = sqrt(x*x + z*z)   inside the plane!
        // |v(sy)| = k(sy) * |v0|
        const int minsy = 40;
        while(lastsy >= minsy) {
            k = scrDistAlongRay / s; // <= 1
            int sy = cast(int) (k * H/2);
            if (sy < lastsy) {
                x = state.data[0]; z = state.data[2];
                int ix = cast(int)x, iz = cast(int)z;
                int c = ((ix ^ iz) & 255) | 128;
                uint clr = (c << 16) + (c << 8) + c;
                img.putPixel(sx + W/2, sy, clr);
                lastsy = sy;
                if (sy==minsy) break;
            }


            double du = state.data[1] * dt, dv = state.data[3] * dt;
            double ds = sqrt(du*du + dv*dv);
            state = evolveRK!geodPlane(state, dt);
            t += dt;
            s += ds;
        }
    }// for sx
}

extern (C) int UIAppMain(string[] args) {
    Window window = Platform.instance.createWindow("Curved", null, 1, 800, 700);
    window.mainWidget = parseML(q{
        VerticalLayout {
            margins: 10pt
            padding: 10pt
            layoutWidth: 750
            HorizontalLayout {
                TextWidget { text: "Hello World" }
                Button {text: "Render"; id: "btnRender"}
            }
            MotionPicture {id: "pic"}
        }
    });

    testRK();
    Vec!(double, 4) v;
    v.data[] = [1.0,2,3,4];
    auto v2 = v * 1.1 + v;
    writeln(v2.data);

    auto imgZ = new ImageZ(512,512);
    auto pic = window.mainWidget.childById!MotionPicture("pic");

    window.mainWidget.childById!Button("btnRender").click = delegate(Widget w) {
        //drawSphere(imgZ);
        drawFloorDirect(imgZ);
        drawFloorRay(imgZ);
        pic.drawImgZ(imgZ);
        return true;
    };
    // you can access loaded items by id - e.g. to assign signal listeners
    // auto edit1 = window.mainWidget.childById!EditLine("edit1");
    // auto edit2 = window.mainWidget.childById!EditLine("edit2");
    // close window on Cancel button click
    // window.mainWidget.childById!Button("btnCancel").click = delegate(Widget w) {
    //     window.close();
    //     return true;
    // };
    // show message box with content of editors
    // window.mainWidget.childById!Button("btnOk").click = delegate(Widget w) {
    //     window.showMessageBox(UIString.fromRaw("Ok button pressed"d),
    //                           UIString.fromRaw("Editors content\nEdit1: "d ~ edit1.text ~ "\nEdit2: "d ~ edit2.text));
    //     return true;
    // };

    // show window
    window.show();

    // run message loop
    return Platform.instance.enterMessageLoop();
}
