//
// Created by edavis on 5/31/19.
//

#pragma once
#ifndef PROTO_PROTO_DATAFLOW_H
#define PROTO_PROTO_DATAFLOW_H

#include <pdfg/GraphIL.hpp>
using namespace pdfg;

using std::cerr;
using std::endl;
using std::map;
using std::ostringstream;
using std::string;
using std::vector;
using std::initializer_list;

// Macros...
#define offset2(i,j,M) ((j)+(i)*(M))
#define offset3(i,j,k,M,N) ((k)+((j)+(i)*(M))*(N))
#define offset4(i,j,k,l,M,N,P) ((l)+((k)+((j)+(i)*(M))*(N))*(P))

namespace Proto {
    struct ConsToPrim : public Macro {
    public:
        explicit ConsToPrim(const Space& output, const Func& func, const vector<Iter>& iters) :
                 Macro(func.name(), iters) {
            Space W = output;
            Space U(func.arg(0).text(), iters);
            Expr gamma = func.arg(1);
            Real half(0.5);
            Real one(1.0);
            Real zero(0.0);

            Iter y = iters[0];
            Iter x = iters[1];
            Iter z;
            if (DIM > 2) {
                z = iters[2];
            }

            // rho = a_U(0);
            // W(0) = rho;
            auto rho = U(Int(0),y,x);
            auto expr = (W(Int(0),y,x) == rho);
            _exprs.push_back(expr);

            // double v2 = 0.0;
            Math v2;
            for (int i = 1; i <= DIM; i++) {
                // double v = a_U(i) / rho;
                Math v = U(i,y,x) / rho;
                // a_W(i) = v;
                expr = (W(i,y,x) == v);
                _exprs.push_back(expr);
                // v2 += v*v;
                auto vsq = paren(paren(v) * paren(v));
                if (v2.empty()) {
                    v2 = vsq;
                } else {
                    v2 = v2 + vsq;
                }
            }
            //a_W(NUMCOMPS-1) = (a_U(NUMCOMPS-1) - half * rho * v2) * (gamma - one);
            expr = (W(Int(NUMCOMPS-1),y,x) == (U(Int(NUMCOMPS-1),y,x) - half * rho * paren(v2)) * paren(gamma - one));
            _exprs.push_back(expr);
        }
    };

    struct WaveSpeedBound : public Macro {
    public:
        explicit WaveSpeedBound(const Space &output, const Func &func, const vector<Iter> &iters) :
                 Macro(func.name(), iters) {
            Space speed = output;
            Space W(func.arg(0).text(), iters);
            Expr gamma = func.arg(1);

            Iter y = iters[0];
            Iter x = iters[1];
            Iter z;
            if (DIM > 2) {
                z = iters[2];
            }

            Math root = Int(DIM) * sqrt(gamma * paren(W(NUMCOMPS-1,y,x) / W(0,y,x)));
            Math expr = (speed(y,x) = root);
            _exprs.push_back(expr);

            Math sum = (speed(y,x) + 0);
            for (int d = 1 ; d <= DIM; d++) {
                //a_speed(0) += a_W(dir);
                sum = sum + W(d,y,x);
            }
            expr = (speed(y,x) = sum);
            _exprs.push_back(expr);
        }
    };

    struct GetFlux : public Macro {
    public:
        explicit GetFlux(const Space &output, const Func &func, const vector<Iter> &iters) :
                Macro(func.name(), iters) {
            Space F = output;
            Space W(func.arg(0).text(), iters);
            Int dim(func.arg(1));
            Real gamma(func.arg(2));

            Iter y = iters[0];
            Iter x = iters[1];
            Iter z;
            if (DIM > 2) {
                z = iters[2];
            }

//            double F0 = a_W(a_d+1)*a_W(0);
//            a_F(0) = F0;
            Math f0 = (F(0,y,x) = W(dim+1,y,x) * W(0,y,x));
            _exprs.push_back(f0);

            Math w2;
            for (int d = 1; d <= DIM; d++) {
//                double Wd = a_W(d);
//                a_F(d) = Wd*F0;
                Math fd = (F(d,y,x) = W(d,y,x) * F(0,y,x));
                _exprs.push_back(fd);
//                W2 += Wd*Wd;
                w2 = w2 + paren(W(d,y,x) * W(d,y,x));
            }
//
//            a_F(a_d+1) += a_W(NUMCOMPS-1);
            Math fdp1 = (F(dim+1,y,x) += W(NUMCOMPS-1,y,x));
            _exprs.push_back(fdp1);
//            a_F(NUMCOMPS-1) = gamma/(gamma - 1.0) * a_W(a_d+1) * a_W(NUMCOMPS-1) + 0.5 * F0 * W2;
            Math fnm1 = gamma/paren(gamma - 1.0) * W(dim+1,y,x) * W(NUMCOMPS-1,y,x) + 0.5 * F(0,y,x) * paren(w2);
            _exprs.push_back(fnm1);
        }
    };

    struct UpwindState : public Macro {
    public:
        explicit UpwindState(const Space &output, const Func &func, const vector<Iter> &iters) :
                Macro(func.name(), iters) {
            Space out = output;
            Space low(func.arg(0).text(), iters);
            Space high(func.arg(1).text(), iters);
            Int dir(func.arg(2));
            Real gamma(func.arg(3));

            Iter y = iters[0];
            Iter x = iters[1];
            Iter z;
            if (DIM > 2) {
                z = iters[2];
            }

            //double rhobar = ((a_low(0) + a_high(0))*.5);
            Math rhobar = paren(low(0,y,x) + high(0,y,x))*.5;
            //double pbar = (a_low(NUMCOMPS-1) + a_high(NUMCOMPS-1))*.5;
            Math pbar = paren(low(NUMCOMPS-1,y,x) + high(NUMCOMPS-1,y,x))*.5;
            //double ubar = (a_low(a_dir+1) + a_high(a_dir+1))*.5;
            Math ubar = paren(low(dir+1,y,x) + high(dir+1,y,x))*.5;
            //double cbar = sqrt(gamma*pbar/rhobar);
            Math cbar = sqrt(gamma*paren(pbar)/paren(rhobar));
            //double pstar = (a_low(NUMCOMPS-1) + a_high(NUMCOMPS-1))*.5 + rhobar*cbar*(a_low(a_dir+1) - a_high(a_dir+1))*.5;
            Math pstar = paren(low(NUMCOMPS-1,y,x) + high(NUMCOMPS-1,y,x))*.5 +paren(rhobar)*
                         paren(cbar)*paren(low(dir+1,y,x) - high(dir+1,y,x))*.5;
            //double ustar = (a_low(a_dir+1) + a_high(a_dir+1))*.5 + (a_low(NUMCOMPS-1) - a_high(NUMCOMPS-1))/(2*rhobar*cbar);
            Math ustar = paren(low(dir+1,y,x) + high(dir+1,y,x))*.5 + paren(low(NUMCOMPS-1,y,x) -
                         high(NUMCOMPS-1,y,x))/paren(2.*paren(rhobar)*paren(cbar));
            Math sign = -sgn(ustar);

            Constr cond = 0. < ustar;
            for (int i = 0; i < NUMCOMPS; i++) {
                Math asgn = (out(i,y,x) = Condition(cond, low(i,y,x), high(i,y,x)) + 0);
                _exprs.push_back(asgn);
            }

            //if (cbar + sign*ubar > 0) {
            Math bar = cbar + sign * ubar;
            cond = 0. < bar;
            //a_out(0) += (pstar - a_out(NUMCOMPS-1))/(cbar*cbar);
            Math diff = paren(pstar - out(NUMCOMPS-1,y,x)) / paren(cbar*cbar);
            Math plus = (out(0,y,x) += Condition(cond, diff, Real(0.)));
            _exprs.push_back(plus);

            //a_out(a_dir+1) = ustar;
            Math outdir = (out(dir+1,y,x) = Condition(cond, ustar, out(dir+1,y,x)) + 0);
            _exprs.push_back(outdir);

            //a_out(NUMCOMPS-1) = pstar;
            Math outnc = (out(NUMCOMPS-1,y,x) = Condition(cond, pstar, out(NUMCOMPS-1,y,x)) + 0);
            _exprs.push_back(outnc);
        }
    };

    class DataFlowFactory {
    public:
        static DataFlowFactory& get() {
            static DataFlowFactory instance;    // Guaranteed to be destroyed.
            return instance;                    // Instantiated on first use.
        }

        DataFlowFactory(DataFlowFactory const&) = delete;         // Don't Implement

        void operator=(DataFlowFactory const&) = delete;     // Don't implement

        void init(const string& name, const string& retname = "", const string& datatype = "",
                  const string& indextype = "") {
            pdfg::init(name);
            if (!datatype.empty()) {
                pdfg::GraphMaker::get().dataType(datatype);
            }
            if (!indextype.empty()) {
                pdfg::GraphMaker::get().indexType(indextype);
            }
            if (!retname.empty()) {
                pdfg::GraphMaker::get().returnName(retname);
            }
        }

        void print(const string& file = "") {
            pdfg::print(file);
        }

        string codegen(const string& path = "") {
            return pdfg::codegen(path);
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        void print(const string& name, const BoxData<T, C, D, E>& boxData) {
            unsigned size = boxData.size();
            const double *data = boxData.data();
            Box box = boxData.box();
            int xmin = (box.low())[0];
            int ymin = (box.low())[1];
            int xmax = (box.high())[0];
            int ymax = (box.high())[1];
            int xsize = (xmax - xmin) + 1;
            int ysize = (ymax - ymin) + 1;

            fprintf(stderr, "%s: [", name.c_str());
            for (unsigned i = 0; i < size; i++) {
                fprintf(stderr, "%lg ", data[i]);
            }
            fprintf(stderr, "]\n");

//            fprintf(stderr, "(c,y,x):\n");
//            for (int c = 0; c < C; c++) {
//                for (int y = ymin; y <= ymax; y++) {
//                    for (int x = xmin; x <= xmax; x++) {
//                        unsigned i = offset3(c, y - ymin, x - xmin, ysize, xsize);
//                        double v = data[i];
//                        fprintf(stderr, "%d,%lg\n", i, v);
//                    }
//                }
//            }
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        Space newSpace(const BoxData<T, C, D, E>& boxData) {
            return newSpace("", boxData);
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        Space newSpace(const string& name, const BoxData<T, C, D, E>& boxData) {
            string sname;
            if (name.empty()) {
                char str[128];
                sprintf(str, "%p", boxData.data());
                sname = string(str);
            } else {
                sname = name;
            }
            _dataPtrs[boxData.data()] = sname;

            vector<Iter> iters;
            vector<pdfg::Range> ranges;
            vector<int> offsets;

            int csizes[3] = {E, D, C};
            char citers[3] = {'e', 'd', 'c'};

            for (unsigned c = 0; c < 3; c++) {
                if (csizes[c] > 1) {
                    Iter iter(citers[c]);
                    iters.emplace_back(iter);
                    pdfg::Range range = (Int(0) <= iter < Int(csizes[c]));
                    ranges.emplace_back(range);
                    offsets.emplace_back(0);
                }
            }

            Box box = boxData.box();
            Point low = box.low();
            Point high = box.high();

            for (int d = 0; d < DIM; d++) {
                Iter iter(_iterNames[d]);
                iters.emplace_back(iter);
                pdfg::Range range = (Int(low[d]) <= iter <= Int(high[d]));
                ranges.emplace_back(range);
                offsets.emplace_back(-low[d]);
            }

            vector<Constr> constrs;
            for (pdfg::Range& range : ranges) {
                constrs.emplace_back(range.lower());
                constrs.emplace_back(range.upper());
            }

            Space space(name, constrs);
            addSpace(space);

            // Create an access for the new space
            vector<Expr> tuple(iters.begin(), iters.end());
            Access access(space, tuple, '(', offsets);

            return space;
        }

        Space newSpace(const Box box) {
            char str[128];
            sprintf(str, "%p", &box);
            return newSpace(string(str), box);
        }

        Space newSpace(const string& name, const Box box) {
            Point low = box.low();
            Point high = box.high();

            vector<Constr> constrs;
            for (int d = 0; d < DIM; d++) {
                Iter iter(_iterNames[d]);
                pdfg::Range range = (Int(low[d]) <= iter <= Int(high[d]));
                constrs.emplace_back(range.lower());
                constrs.emplace_back(range.upper());
            }

            Space space(name, constrs);
            addSpace(space);
            return space;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Macro newMacro(const string& name, Space& input, Space& output, const pdfg::Func& func) {
            Macro macro;
            if (name.find("consToPrim") != string::npos) {
                macro = ConsToPrim(output, func, input.iterators());
            } else if (name.find("waveSpeedBound") != string::npos) {
                macro = WaveSpeedBound(output, func, input.iterators());
            } else if (name.find("getFlux") != string::npos) {
                macro = GetFlux(output, func, input.iterators());
            } else if (name.find("upwindState") != string::npos) {
                macro = UpwindState(output, func, input.iterators());
            } else {
                int stop = 1;
            }
            return macro;
        }

//        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1, typename Func, typename... Srcs>
//        Comp newComp(const Func& func, Srcs&... srcs) {
//            ostringstream os;
//            os << &func;
//            return newComp(os.str(), func, srcs...);
//        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1, typename Func, typename... Srcs>
        Comp newComp(const string& name, initializer_list<string> inNames, const string& outName,
                     const BoxData<T, C, D, E>& outData, const Func& func, Srcs&... srcs) {
            Box box;
            return newComp(name, inNames, outName, outData, func, box, srcs...);
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1, typename Func, typename... Srcs>
        Comp newComp(const string& name, initializer_list<string> inNames, const string& outName,
                     const BoxData<T, C, D, E>& outData, const Func& func, Box domain, Srcs&... srcs) {
            // Copy input names...
            _inputPos = 0;
            _inputNames.clear();
            for (const string& input : inNames) {
                _inputNames.emplace_back(input);
            }

            // Compute intersection to get iteration space for 'forall'.
            Box intersect;
            pdfg::Func fxn(name);
            parseArgs(intersect, fxn, srcs...);

            if (domain.empty()) {
                domain = intersect;
            }

            // Create input/output spaces.
            Space ispace = newSpace("I" + name, domain);
            Space ospace = newSpace(outName, outData);

            Macro macro = newMacro<T,C,D,E>(name, ispace, ospace, fxn);
            Comp comp(macro.name(), ispace, call(macro));
            return comp;
        }

//        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
//        Comp newComp(const string& name, const string& outName, const Stencil<T>& stencil,
//                     const BoxData<T, C, D, E>& inData, const BoxData<T, C, D, E>& outData) {
//            return (name, "", outName, stencil, inData, outData);
//        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const Stencil<T>& stencil,
                     const BoxData<T, C, D, E>& inData, T scale, const BoxData<T, C, D, E>& outData) {
            auto itr = _dataPtrs.find(inData.data());
            Space dspace;
            if (itr != _dataPtrs.end()) {
                dspace = getSpace(itr->second);
            } else {
                dspace = newSpace("", inData);
            }

            Box domain = inData.box();
            Space ispace = newSpace("I" + name, domain);
            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            // TODO: Calculate output space size for stencil...
            Space ospace = newSpace(outName, outData);
            Access lhs(ospace, tuple);

            vector<T> coefs = stencil.coefs();
            vector<Point> points = stencil.offsets();
            Math rhs;
            for (unsigned k = 0; k < coefs.size(); k++) {
                Real coef(coefs[k]);
                Point point = points[k];

                tuple.clear();
                for (int d = 0; d < DIM; d++) {
                    tuple.push_back(iters[d] + point[d]);
                }

                Access acc(dspace, tuple);
                Math prod = coef * acc;
                if (k < 1) {
                    rhs = prod;
                } else {
                    rhs = rhs + prod;
                }
            }

            // TODO: Figure out how to get scale factor, to multiply result by...
//            if (scale == 1.0) {
//                scale /= (T) coefs.size();
//            }

            Math expr = lhs = rhs * scale;
            Comp comp(name, ispace, (expr));
            return comp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const string& oper,
                     T lhs, const BoxData<T, C, D, E>& rhs) {
            Box domain = rhs.box();
            Space lspace(outName);
            Space rspace = getSpace(_dataPtrs[rhs.data()]);

            Space ispace = newSpace("I" + name, domain);
            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            Access readAcc(rspace, tuple);
            Math expr(lspace, readAcc, oper);

            // TODO: This results in a cyclic reference, need to resolve
//            if (oper.find('(') != string::npos) {
//                Math assign(writeAcc, expr, "=");       // Assign write access to result of function...
//                expr = assign;
//            }

            Comp comp(name, ispace, (expr));
            return comp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const string& oper,
                     const BoxData<T, C, D, E>& lhs, T rhs) {
            Box domain = lhs.box();
            Space lspace = getSpace(_dataPtrs[lhs.data()]);

            Space ispace = newSpace("I" + name, domain);
            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            Access writeAcc(lspace, tuple);
            Real scalar(rhs);
            Math expr(writeAcc, scalar, oper);

            Comp comp(name, ispace, (expr));
            return comp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const string& oper,
                     const BoxData<T, C, D, E>& lhs, const BoxData<T, C, D, E>& rhs) {
            Box lbox = lhs.box();
            Box rbox = rhs.box();
            Box domain;
            if (lbox == rbox) {
                domain = lbox;
            } else {
                domain = lbox & rbox;
            }

            //Space ospace = newSpace(outName, lhs);
            Space lspace = getSpace(_dataPtrs[lhs.data()]);
            Space rspace = getSpace(_dataPtrs[rhs.data()]);

            return newComp<T,C,D,E>(name, outName, oper, domain, lspace, rspace);
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& opName, const string& outName, const string& oper,
                     const BoxData<T, C, D, E>& lhs, const BoxData<T, C, D, E>& rhs,
                     const string& stencilName, const Stencil<T>& stencil, T scale) {
            Box lbox = lhs.box();
            Box rbox = rhs.box();
            Box domain;
            if (lbox == rbox) {
                domain = lbox;
            } else {
                domain = lbox & rbox;
            }

            // 1) Apply stencil computation with temporary intermediate output.
//            Comp newComp(const string& name, const string& outName, const Stencil<T>& stencil,
//                         const BoxData<T, C, D, E>& inData, T scale, const BoxData<T, C, D, E>& outData)
            string stencilOut = stencilName + "_out";
            Space tspace = newSpace(stencilOut, lhs);
            Comp stencilComp = newComp(stencilName, stencilOut, stencil, rhs, scale, lhs);

            // 2) Apply operation...
            //Space lspace = getSpace(_dataPtrs[lhs.data()]);
            Space rspace = getSpace(_dataPtrs[rhs.data()]);

            return newComp<T,C,D,E>(opName, outName, oper, domain, tspace, rspace);
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const string& oper,
                     Box domain, Space lspace, Space rspace) {
            Space ispace = newSpace("I" + name, domain);
            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            Access writeAcc(lspace, tuple);
            Access readAcc(rspace, tuple);
            Math expr(writeAcc, readAcc, oper);

            Comp comp(name, ispace, (expr));
            return comp;
        }

    protected:
        DataFlowFactory() = default;

        template <typename T, unsigned int C, unsigned char D, unsigned char E>
        void parse(Box& domain, Func& func, const BoxData<T,C,D,E>& data) {
            // Recursively compute domain intersection...
            if (domain.empty()) {
                domain = data.box();
            } else {
                domain = domain & data.box();
            }

            // Find or create corresponding data space.
            Space space;
            auto itr =  _dataPtrs.find(data.data());
            if (itr != _dataPtrs.end()) {
                space = getSpace(itr->second);
            } else {
                space = newSpace(_inputNames[_inputPos], data);
                _inputPos += 1;
            }

            // Add data space to function as arg...
            func.add(space);
        }

        template <typename T>
        void parse(Box& domain, Func& func, T scalar) {
            // Satisfy compiler, we can convert this into a 'Real' object, but then what?
            Real real(scalar);
            func.add(real);
            //_scalars.push_back(real);
        }

        template <typename Last>
        void parseArgs(Box& domain, Func& func, Last& last) {
            parse(domain, func, last);
        }

        template <typename T, unsigned int C, unsigned char D, unsigned char E, typename... Rest>
        void parseArgs(Box& domain, Func& func, const BoxData<T,C,D,E>& first, Rest&... rest) {
            parse(domain, func, first);
            parseArgs(domain, func, rest...);
        }

        template <typename T, unsigned int C, unsigned char D, unsigned char E>
        void parseArgs(Box& domain, Func& func, const BoxData<T,C,D,E>& first) {
            parse(domain, func, first);
        }

        template <typename First , typename... Rest>
        void parseArgs(Box& domain, Func& func, First& first, Rest&... rest) {
            parse(domain, func, first);
            parseArgs(domain, func, rest...);
        }

#if DIM>2
        const char _iterNames[9] = {'z', 'y', 'x', 'i', 'j', 'k', 'l', 'm', 'n'};
#elif DIM>1
        const char _iterNames[8] = {'y', 'x', 'i', 'j', 'k', 'l', 'm', 'n'};
#else
        const char _iterNames[7] = {'x', 'i', 'j', 'k', 'l', 'm', 'n'};
#endif
        map<const void*, string> _dataPtrs;
        vector<string> _inputNames;
        size_t _inputPos;
    };
}

#endif //PROTO_PROTO_DATAFLOW_H
