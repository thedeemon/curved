import std.stdio, dlangui, dlangui.widgets.metadata, std.math, std.conv, std.array, std.format, std.meta;
import symbolic, img, geometry;
import std.algorithm : canFind;

mixin APP_ENTRY_POINT;
enum VX = 854 - 320;

class DrawingBoard : ImageWidget {
	ColorDrawBuf cdbuf;
    int W = 854, H = 480;

	this() {
		cdbuf = new ColorDrawBuf(W, H);
		Ref!DrawBuf r = cdbuf;
		this.drawable = new ImageDrawable(r);
		this.margins = 0;// Rect(0, 10, 0, 10);//10
        fillBytes(0);
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

        enum mapData = cast(immutable ubyte[]) import("map.png");
        worldMap = loadImage(mapData, "map.png");
    }

    version(DirectorsCut) {
        void save(string fname) { saveBuf(fname, cdbuf); }
        void delegate() timerF;
        override bool onTimer(ulong id) {
            if (timerF !is null)
                timerF();
            return true;
        }
    }
}

mixin(registerWidgets!("__gshared static this", DrawingBoard));
double mod360(double x) { while (x>=360) x -= 360; return x; }

extern (C) int UIAppMain(string[] args) {
    import dlangui.core.logger;
    int w= 880, h = 600;

   	Log.setLogLevel( dlangui.core.logger.LogLevel.Error );

    // foreach(s; genCode!(Wormhole.equation)) writeln(s);
    // return 0;

    version(Windows) {
        w = w.pixelsToPoints; h = h.pixelsToPoints;
    }
    Window window = Platform.instance.createWindow("Curved spaces", null, 1, w, h);
    window.mainWidget = parseML(q{
        VerticalLayout {
            backgroundColor: 0xc0c0c0
            padding: 10pt
            layoutWidth: 750
            HorizontalLayout {
                Button {text: "Render"; id: "btnRender"}
                ComboBox {id: "world"}
                TextWidget {text: "Heading:" }
                EditLine { text: "90"; id: "heading"; layoutWidth: 50}
                // TextWidget {text: "dt:" }
                // EditLine { text: "1"; id: "dt"; layoutWidth: 50}
                CheckBox { text: "Dyn. step"; id: "dyndt"}
                CheckBox { text: "Walls"; id: "walls"}
                TextWidget {text: "Range:" }
                EditLine { text: "8"; id: "range"; layoutWidth: 70}
                TextWidget {text: "R:" }
                EditLine { text: "1000"; id: "R"; layoutWidth: 70}
                TextWidget {text:""; id:"out"}
            }
            DrawingBoard {id: "pic"}
        }
    });

    auto edHeading = window.mainWidget.childById!EditLine("heading");
    auto txtOut = window.mainWidget.childById!TextWidget("out");
    auto cbDDT = window.mainWidget.childById!CheckBox("dyndt");
    auto cbWalls = window.mainWidget.childById!CheckBox("walls");
    auto edRange = window.mainWidget.childById!EditLine("range");
    auto edR = window.mainWidget.childById!EditLine("R");
    auto worldCombo = window.mainWidget.childById!ComboBox("world");
    auto pic = window.mainWidget.childById!DrawingBoard("pic");
    auto vl = cast(VerticalLayout) window.mainWidget;

    version(DirectorsCut) {
        int shotN = 0, batchFramesMove = 0, batchFramesView = 0;
        auto ui = parseML(q{
            HorizontalLayout {
                EditLine { text: "/home/dee/tmp/curved/one/"; id: "picdir"; layoutWidth: 400 }
                Button { text: "Save"; id: "saveBtn" }
                CheckBox { text: "action!"; id: "record" }
                TextWidget { text: "..."; id: "times" }
            }
        });
        pic.setTimer(10);
        auto edPath = ui.childById!EditLine("picdir");
        auto cbRecord = ui.childById!CheckBox("record");
        auto txtTimes = ui.childById!TextWidget("times");
        vl.addChild(ui);

        void saveFrame() {
            import std.file : exists, mkdir;
            import std.path : buildPath;
            auto path = edPath.text;
            if (!exists(path))
                mkdir(path);
            pic.save( buildPath(path.to!string, format("%05d.png", shotN++)) );
        }

        ui.childById!Button("saveBtn").click = delegate(Widget w) {
            saveFrame();
            return true;
        };
    } else {
        vl.addChild(new TextWidget().text("WASD - move and turn, IJKL - rotate 3D view, U/O - zoom it."));
    }

    Renderer[] worlds;
    foreach(Wld; AliasSeq!(Earth, Sphere, FatBall, Donut, BlackHole, Wormhole))
        worlds ~= new Render!Wld;

    worldCombo.items = worlds.amap!(w => w.name.to!dstring);
    auto imgSphere = new ImageZ(512,512);
    auto img = new ImageZ(512,512);
    auto imgFloor = new ImageZ(VX, 480);

    auto ps = new Params();
    ps.pos.u = 4.7; ps.pos.v = 0.1;
    ps.range = 8.0; ps.dt = 1;
    ps.rotAlpha = 0; ps.rotBeta = 0;
    Renderer rend = worlds[0];
    rend.setDistance();
    rend.drawSurface(imgSphere, ps);
    UV[] points;

    void render() {
        import std.datetime.stopwatch : StopWatch, AutoStart;
        auto sw = StopWatch(AutoStart.yes);
        ps.heading = edHeading.text.to!double;
        ps.dyndt = cbDDT.checked;
        ps.walls = cbWalls.checked;
        auto rng = edRange.text.to!double;
        if (rng >= 1 && rng <= 20) ps.range = rng;
        auto r = edR.text.to!double;
        if (r >= 100 && r <= 1000000) globalR = r;
        img.copyFrom(imgSphere);
        points = rend.drawFloorRay(imgFloor, ps, ps.pos.u, ps.pos.v);
        rend.drawPoints(img, points, ps);
        pic.drawImgAt(img, 0,0, 100,100, 320,320);
        pic.drawImgAt(imgFloor, 320,0, 0,0, VX,480);
        version(DirectorsCut)
            txtOut.text = format("t=%s batchM=%s batchV=%s"d, sw.peek.total!"msecs",
                                    batchFramesMove, batchFramesView);
        else
            txtOut.text = format("t=%s"d, sw.peek.total!"msecs");
    }

    void rerender3DView() {
        rend.drawSurface(imgSphere, ps);
        img.copyFrom(imgSphere);
        rend.drawPoints(img, points, ps);
        pic.drawImgAt(img, 0,0, 100,100, 320,320);
        window.invalidate();
    }

    window.mainWidget.childById!Button("btnRender").click = delegate(Widget w) {
        render();
        return true;
    };

    worldCombo.itemClick = delegate(Widget wgt, int idx) {
        rend = worlds[idx];
        rend.setDistance();
        rend.ensureCorrectPos(ps);
        rend.drawSurface(imgSphere, ps);
        render();
        return true;
    };

    bool moveAndDraw(uint key, double k) {
        switch(key) {
            case KeyCode.KEY_A: ps.heading -= 10*k; edHeading.text = ps.heading.to!dstring; break;
            case KeyCode.KEY_D: ps.heading += 10*k; edHeading.text = ps.heading.to!dstring; break;
            case KeyCode.KEY_W:
                auto hdn = ps.heading;
                ps.pos = rend.walk(ps.pos, hdn, ps.dt, 50 * k);
                ps.heading = hdn;
                edHeading.text = format("%.3g"d, hdn);
                break;
            case KeyCode.KEY_S:
                auto hdn = ps.heading + 180;
                ps.pos = rend.walk(ps.pos, hdn, ps.dt, 50 * k);
                ps.heading = hdn - 180;
                while(ps.heading < 0) ps.heading += 360;
                while(ps.heading >= 360) ps.heading -= 360;
                edHeading.text = format("%.3g"d, ps.heading);
                break;
            case KeyCode.KEY_L: ps.rotAlpha = mod360(ps.rotAlpha + 360-10*k); rerender3DView(); return true;
            case KeyCode.KEY_J: ps.rotAlpha = mod360(ps.rotAlpha + 10*k);     rerender3DView(); return true;
            case KeyCode.KEY_I: ps.rotBeta = mod360(ps.rotBeta + 360-10*k);   rerender3DView(); return true;
            case KeyCode.KEY_K: ps.rotBeta = mod360(ps.rotBeta + 10*k);       rerender3DView(); return true;
            case KeyCode.KEY_U: globalDistance += 0.5*k; rerender3DView(); return true;
            case KeyCode.KEY_O: globalDistance -= 0.5*k; rerender3DView(); return true;
            default: return false;
        }
        render();
        window.invalidate();
        return true;
    }

    version(DirectorsCut) {
        uint last_key_move, last_key_view;

        void film() {
            if (batchFramesMove + batchFramesView <= 0) return;
            if (batchFramesMove > 0) {
                moveAndDraw(last_key_move, 0.03);
                batchFramesMove--;
            }
            if (batchFramesView > 0) {
                moveAndDraw(last_key_view, 0.05);
                batchFramesView--;
            }
            saveFrame();
        }

        pic.timerF = { if (cbRecord.checked) film(); };
    }

    window.mainWidget.keyEvent = delegate(Widget wt, KeyEvent e) {
        //writeln("main keyEvent ", e);
        if (e.action==KeyAction.KeyDown) {
            version(DirectorsCut) {
                static keysMove = [KeyCode.KEY_W, KeyCode.KEY_A, KeyCode.KEY_S, KeyCode.KEY_D];
                static keysView = [KeyCode.KEY_U, KeyCode.KEY_I, KeyCode.KEY_O, KeyCode.KEY_J, KeyCode.KEY_K, KeyCode.KEY_L];
                if (cbRecord.checked) {
                    if (keysMove.canFind(e.keyCode)) {
                        batchFramesMove += 30;
                        last_key_move = e.keyCode;
                    }
                    if (keysView.canFind(e.keyCode)) {
                        batchFramesView += 20;
                        last_key_view = e.keyCode;
                    }
                    return true;
                }
            }
            return moveAndDraw(e.keyCode, 1);
        }
        return false;
    };

    window.show();
    render();
    return Platform.instance.enterMessageLoop();
}
