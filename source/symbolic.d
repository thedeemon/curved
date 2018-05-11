import std.stdio, std.conv, std.format, std.array : array, replace;
import std.algorithm : map, filter, sort, all, any, canFind, count, countUntil;
import std.string : join;
import std.range : iota;

class Expr {
    abstract Expr diff(string dx);
    abstract string code();
}

class Const : Expr {
    string v; // must contain a number
    this(string x) { v = x; }
    override string toString() {  return v;  }
    override string code() { return v[0]=='-' ? "("~v~")" : v; }
    override Expr diff(string dx) { return new Const("0"); }
}

class Var : Expr {
    string var;
    this(string v) { var = v; }
    override string toString() { return var; }
    override string code() { return var; }
    override Expr diff(string dx) {
        if (dx==var) return new Const("1");
        return zero;
    }
}

class Sin : Expr { 
    Expr arg;
    this(string v) { arg = new Var(v); }
    this(Expr x) { arg = x; }
    override string toString() { return format("sin(%s)", arg); }
    override string code() { return "sin_"~ arg.code; }
    override Expr diff(string dx) {
        return mul(new Cos(arg), arg.diff(dx));
    }
}

class Cos : Expr {
    Expr arg;
    this(string v) { arg = new Var(v); }
    this(Expr x) { arg = x; }
    override string toString() { return format("cos(%s)", arg); }
    override string code() { return "cos_"~ arg.code; }
    override Expr diff(string dx) {
        return mul(neg(new Sin(arg)), arg.diff(dx));
    }
}

class Mul : Expr {
    Expr[] abc;
    this(Expr x, Expr y) { abc =[x,y]; }
    this(Expr[] xs) {
        foreach(x; xs) {
            if (Mul m = cast(Mul)x) abc ~= m.abc;
            else abc ~= x;
        }
    }

    override string toString() { return "["~ abc.map!str.join(" * ") ~"]"; }
    override string code() { return abc.map!toCode.join(" * "); }
    override Expr diff(string dx) {
        Expr[] added;
        foreach(i; 0..abc.length) {
            auto di = abc[i].diff(dx);
            Expr[] xs;
            if (Mul m = cast(Mul)di)
                xs = abc[0..i] ~ m.abc ~ abc[i+1..$];
            else
                xs = abc[0..i] ~ di ~ abc[i+1..$];
            added ~= new Mul(xs).simp;
        }
        return new Add(added).simp;
    }

    Expr simp() {
        if (abc.canFind!(e => e.eq(0))) return zero;
        sort!((a,b) => a.toString < b.toString)(abc);
        auto nonTriv = abc.filter!(e => !e.eq(1)).array;
        if (nonTriv.length >= 2 && nonTriv[0].eq(-1) && nonTriv[1].eq(-1))
            nonTriv = nonTriv[2..$];
        Expr[] nums, denums, nonDivs;
        foreach(e; nonTriv) {
            if (Div d = cast(Div)e) {
                nums ~= d.a; denums ~= d.b;
            } else
                nonDivs ~= e;
        }
        if (nums.length > 0)
            return new Div( new Mul(nonDivs ~ nums).simp, new Mul(denums).simp).simp;

        if (nonTriv.length==abc.length) return this;
        switch(nonTriv.length) {
            case 0: return new Const("1");
            case 1: return nonTriv[0];
            default: return new Mul(nonTriv);
        }
    }
}

class Add : Expr {
    Expr[] abc;
    this(Expr x, Expr y) { abc = [x,y]; }
    this(Expr[] xs) { abc = xs; }
    override string toString() { return "("~ abc.map!str.join(" + ") ~")"; }
    override string code() { return "(" ~ abc.map!toCode.join(" + ") ~")"; }
    override Expr diff(string dx) {
        return new Add( abc.map!(e => e.diff(dx)).array ).simp;
    }

    Expr simp() {
        sort!((a,b) => a.toString < b.toString)(abc);
        auto nonZeroes = abc.filter!(e => !e.eq(0)).array;

        if (nonZeroes.length > 1) { // try factoring
            nonZeroes = pythagoras(nonZeroes);

            while(nonZeroes.length > 1) {
                Const c0 = cast(Const)nonZeroes[0], c1 = cast(Const)nonZeroes[1];
                if (c0 && c1) {
                    auto c = c0.v.to!int + c1.v.to!int;
                    Expr e = new Const(c.to!string);
                    nonZeroes = [e] ~ nonZeroes[2..$];
                } else break;
            }
            if (nonZeroes.length > 1) {
                bool ok = false;
                Expr factored = tryFactoring(nonZeroes, ok);
                if (ok) return factored;
            }
        }

        if (nonZeroes.length == abc.length) return this; // nothing changed
        switch(nonZeroes.length) {
            case 0: return zero;
            case 1: return nonZeroes[0];
            default: return new Add( nonZeroes );
        }
    }
}

class Div : Expr {
    Expr a, b; // a/b
    this(Expr x, Expr y) { a=x; b=y; }
    override string toString() { return "{" ~ a.str ~ " / " ~ b.str ~ "}"; }
    override string code() { return "(("~a.code~") / ("~b.code~"))"; }
    override Expr diff(string dx) {
        auto da = a.diff(dx), db = b.diff(dx);
        return  new Div( add( mul(da,b), neg(mul(db, a)) ),  mul(b,b)).simp;
    }

    Expr simp() {
        Expr[] fa = factors(a), fb = factors(b);
        Expr[] common = commonFactors(fa, fb);
        if (common.length > 0) {
            string[] coms = common.amap!str;
            auto na = factorOut(fa, coms), nb = factorOut(fb, coms);
            return divSimp(na, nb);
        }
        return divSimp(a,b);
    }
}

Expr binOp(F)(Expr a, Expr b) { //  Mul or Add, flatten and simplify
    F ma = cast(F)a, mb = cast(F)b;
    if (ma && mb) return new F(ma.abc ~ mb.abc).simp;
    if (ma) return new F(ma.abc ~ b).simp;
    if (mb) return new F(mb.abc ~ a).simp;
    return new F(a,b).simp;
}

Expr mul(Expr a, Expr b) { return binOp!Mul(a,b); }
Expr add(Expr a, Expr b) { return binOp!Add(a,b); }
Expr mul(Expr[] xs...) { return new Mul(xs).simp; }
Expr add(Expr[] xs...) { return new Add(xs).simp; }
Expr div(Expr a, Expr b) {
    Div da = cast(Div)a, db = cast(Div)b;
    if (da) {
        if (db) return div(mul(da.a, db.b), mul(da.b, db.a));
        else return new Div(da.a, mul(da.b, b)).simp;
    } else {
        if (db) return new Div(mul(a,db.b), db.a).simp;
    }
    return new Div(a,b).simp;
}

auto amap(alias f, R)(R xs) { return xs.map!f.array; }
Const zero() { return new Const("0"); }
Expr neg(Expr a) { return mul(new Const("-1"), a); }
string str(Expr a) { return a.toString; }
string toCode(Expr a) { return a.code; }
bool eq(Expr a, Expr b) { return a.toString == b.toString; }
bool eq(Expr a, string s) { return a.toString == s; }
bool eq(Expr e, int x) {
    if (Const c = cast(Const)e)
        return c.v==x.to!string;
    return false;
}

Expr[] factors(Expr e) {
    if (Mul m = cast(Mul)e) {
        Expr[] arr = m.abc;
        sort!((a,b) => a.toString < b.toString)(arr);
        return arr;
    }
    return [e];
}

string[] factorNames(Expr e) {
    return e.factors.amap!str;
}

Expr[] commonFactors(Expr[] xs, Expr[] ys) {
    Expr[] common;
    size_t i=0, j=0;
    while(i < xs.length && j < ys.length) {
        if (xs[i].eq(ys[j])) {
            common ~= xs[i];
            i++; j++;
        } else {
            if (xs[i].toString < ys[j].toString) i++;
            else j++;
        }
    }
    return common;
}

Expr factorOut(Expr[] fs, string[] common) {
    Expr[] arr;
    size_t i=0, j=0;
    while(i < fs.length && j < common.length) {
        if (fs[i].eq(common[j])) {
            i++; j++;
        } else {
            if (fs[i].toString < common[j]) {
                arr ~= fs[i];
                i++;
            }
            else j++;
        }
    }
    if (i < fs.length) arr ~= fs[i..$];

    switch(arr.length) {
        case 0: return new Const("1");
        case 1: return arr[0];
        default: return new Mul(arr).simp;
    }
}

bool divisibleBy(string[] allFactors, string[] by) {
    size_t i=0, j=0;
    while(i < allFactors.length && j < by.length) {
        if (allFactors[i] == by[j]) {
            i++; j++;
        } else {
            if (allFactors[i] < by[j]) i++;
            else return false;//j++;
        }
    }
    return j==by.length; //true if we exhausted the list
}

struct Divisor {
    string[] names;
    Expr[] exprs;
}

Expr tryFactoring(Expr[] summands, ref bool ok) {
    Expr[][] fs = summands.amap!factors;
    string[][] ns = fs.amap!(a => a.amap!str);
    Divisor[] candidates;
    foreach(i; 0..summands.length-1)
        foreach(j; i+1 .. summands.length) {
            Expr[] a = fs[i]; string[] bs = ns[j];
            Expr[] common = commonFactors(fs[i], fs[j]);
            if (common.length > 0) {
                string[] coms = common.amap!str;
                if (!candidates.canFind!(c => c.names == coms))
                    candidates ~= Divisor(coms,  common);
            }
        }

    if (candidates.length == 0) {
        ok = false; 
        return new Add(summands); // nothing found
    }
    ok = true;
    string[] bestCand;
    size_t bestCount = 0;
    Expr[] bestDivisor;
    foreach(div; candidates) {
        auto n = ns.count!(f => f.divisibleBy(div.names));
        if (n > bestCount) {
            bestCount = n; bestCand = div.names; bestDivisor = div.exprs;
        }
    }

    Expr[] divided, undivided;
    foreach(i, e; summands) {
        if (ns[i].divisibleBy(bestCand)) {
            Expr rest = factorOut(fs[i], bestCand);
            divided ~= rest;
        } else undivided ~= e;
    }
    auto factored = new Mul(bestDivisor ~ (new Add(divided).simp)).simp;

    if (undivided.length > 0)
        return new Add(undivided ~ factored).simp;
    else
        return factored;
}

bool isTrig2(F)(Expr e, string var) { // F is Cos or Sin
    if (Mul m = cast(Mul)e) {
        if (m.abc.length==2) {
            auto c1 = cast(F)m.abc[0], c2 = cast(F)m.abc[1];
            return (c1 && c2 && c1.arg.eq(var) && c2.arg.eq(var));
        }
    }
    return false;
}

bool isCos2(Expr e, string var) { return isTrig2!Cos(e, var); }
bool isSin2(Expr e, string var) { return isTrig2!Sin(e, var); }

// cos(x)^2 + sin(x)^2 = 1
Expr[] pythagoras(Expr[] summands) {
    foreach(var; ["u","v"]) {
        auto cos2 = summands.countUntil!(e => e.isCos2(var));
        auto sin2 = summands.countUntil!(e => e.isSin2(var));
        if (cos2 >= 0 && sin2 >= 0) {
            Expr[] res = [new Const("1")];
            foreach(i; 0..summands.length)
                if (i != sin2 && i != cos2)
                    res ~= summands[i];
            summands = res;
        }
    }
    return summands;
}

Expr divSimp(Expr a, Expr b) {
    if (b.eq(1)) return a;
    if (a.eq(0)) return zero;
    return new Div(a, b);
}

Expr[] diff(Expr[] vec, string dx) { return vec.amap!(e => e.diff(dx)); }
Expr dot(Expr[] a, Expr[] b) { return iota(a.length).amap!(i => mul(a[i], b[i])).add; }

string txtSimp(string s) {
    return s.replace("sqrt(R * R)","R").replace("sqrt(1)", "1");
}