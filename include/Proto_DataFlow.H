//
// Created by edavis on 5/31/19.
//

#pragma once
#ifndef PROTO_PROTO_DATAFLOW_H
#define PROTO_PROTO_DATAFLOW_H

#include <pdfg/GraphIL.hpp>
using namespace pdfg;

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
            vector<Expr>& tuple = _args;

            // rho = a_U(0);
            // W(0) = rho;
            Access rho = U(0,tuple);
            auto expr = (W(0,tuple) == rho);
            _exprs.push_back(expr);

            // double v2 = 0.0;
            Math v2;
            for (int d = 1; d <= DIM; d++) {
                // double v = a_U(d) / rho;
                Math v = U(d,tuple) / rho;
                // a_W(d) = v;
                expr = (W(d,tuple) == v);
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
            expr = (W(NUMCOMPS-1,tuple) == (U(NUMCOMPS-1,tuple) - half * rho * paren(v2)) * paren(gamma - one));
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
            vector<Expr>& tuple = _args;

            //a_speed(0) = DIM*sqrt(gamma*a_W(NUMCOMPS-1)/a_W(0));
            Math sum = paren(Real(DIM) * sqrt(gamma * W(NUMCOMPS-1,tuple) / W(0,tuple)));
            for (int d = 1 ; d <= DIM; d++) {
                //a_speed(0) += a_W(dir);
                sum = sum + W(d,tuple);
            }

            Math asgn = (speed(tuple) = sum);
            _exprs.push_back(asgn);
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
            vector<Expr>& tuple = _args;

//            double F0 = a_W(a_d+1)*a_W(0);
//            a_F(0) = F0;
            Math f0 = (F(0,tuple) = W(dim+1,tuple) * W(0,tuple));
            _exprs.push_back(f0);

            Math w2;
            for (int d = 1; d <= DIM; d++) {
//                double Wd = a_W(d);
//                a_F(d) = Wd*F0;
                Math fd = (F(d,tuple) = W(d,tuple) * F(0,tuple));
                _exprs.push_back(fd);
//                W2 += Wd*Wd;
                w2 = w2 + paren(W(d,tuple) * W(d,tuple));
            }
//
//            a_F(a_d+1) += a_W(NUMCOMPS-1);
            Math fdp1 = (F(dim+1,tuple) += W(NUMCOMPS-1,tuple));
            _exprs.push_back(fdp1);
//            a_F(NUMCOMPS-1) = gamma/(gamma - 1.0) * a_W(a_d+1) * a_W(NUMCOMPS-1) + 0.5 * F0 * W2;
            Math fnm1 = (F(NUMCOMPS-1,tuple) = gamma/paren(gamma - 1.0) * W(dim+1,tuple) * W(NUMCOMPS-1,tuple) +
                                               0.5 * F(0,tuple) * paren(w2));
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
            vector<Expr>& tuple = _args;

            //double rhobar = ((a_low(0) + a_high(0))*.5);
            Math rhobar = paren(low(0,tuple) + high(0,tuple))*.5;
            //double pbar = (a_low(NUMCOMPS-1) + a_high(NUMCOMPS-1))*.5;
            Math pbar = paren(low(NUMCOMPS-1,tuple) + high(NUMCOMPS-1,tuple))*.5;
            //double ubar = (a_low(a_dir+1) + a_high(a_dir+1))*.5;
            Math ubar = paren(low(dir+1,tuple) + high(dir+1,tuple))*.5;
            //double cbar = sqrt(gamma*pbar/rhobar);
            Math cbar = sqrt(gamma*paren(pbar)/paren(rhobar));
            //double pstar = (a_low(NUMCOMPS-1) + a_high(NUMCOMPS-1))*.5 + rhobar*cbar*(a_low(a_dir+1) - a_high(a_dir+1))*.5;
            Math pstar = paren(low(NUMCOMPS-1,tuple) + high(NUMCOMPS-1,tuple))*.5 +paren(rhobar)*
                         paren(cbar)*paren(low(dir+1,tuple) - high(dir+1,tuple))*.5;
            //double ustar = (a_low(a_dir+1) + a_high(a_dir+1))*.5 + (a_low(NUMCOMPS-1) - a_high(NUMCOMPS-1))/(2*rhobar*cbar);
            Math ustar = paren(low(dir+1,tuple) + high(dir+1,tuple))*.5 + paren(low(NUMCOMPS-1,tuple) -
                         high(NUMCOMPS-1,tuple))/paren(2.*paren(rhobar)*paren(cbar));
            Math sign = -sgn(ustar);

            Constr cond = 0. < ustar;
            for (int i = 0; i < NUMCOMPS; i++) {
                Math asgn = (out(i,tuple) = Condition(cond, low(i,tuple), high(i,tuple)) + 0);
                _exprs.push_back(asgn);
            }

            //if (cbar + sign*ubar > 0) {
            Math bar = cbar + sign * ubar;
            cond = 0. < bar;
            //a_out(0) += (pstar - a_out(NUMCOMPS-1))/(cbar*cbar);
            Math diff = paren(pstar - out(NUMCOMPS-1,tuple)) / paren(cbar*cbar);
            Math plus = (out(0,tuple) += Condition(cond, diff, Real(0.)));
            _exprs.push_back(plus);

            //a_out(a_dir+1) = ustar;
            Math outdir = (out(dir+1,tuple) = Condition(cond, ustar, out(dir+1,tuple)) + 0);
            _exprs.push_back(outdir);

            //a_out(NUMCOMPS-1) = pstar;
            Math outnc = (out(NUMCOMPS-1,tuple) = Condition(cond, pstar, out(NUMCOMPS-1,tuple)) + 0);
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
                  const string& indextype = "", const string& defval = "") {
            pdfg::init(name, retname, datatype, indextype, {}, defval);
            pdfg::align_iters();
        }

        void constants(initializer_list<string> names, initializer_list<int> values = {}) {
            pdfg::addConstants(names, values);
        }

        void print(const string& file = "") {
            pdfg::print(file);
        }

        void tile(unsigned size = 8) {
            pdfg::addConstant("T", size);
#if DIM<3
            pdfg::tile({"y", "x"}, {size, size});
#else
            pdfg::tile({"z", "y", "x"}, {size, size, size});
#endif
        }

        string codegen(const string& path = "") {
            pdfg::data_reduce();
            pdfg::mem_alloc();
            return pdfg::codegen(path, "", "C++");
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        void print(const string& name, const BoxData<T, C, D, E>& boxData, const string file = "") {
            unsigned size = boxData.size();
            const double *data = boxData.data();
            Box box = boxData.box();
            int xmin = (box.low())[0];
            int ymin = (box.low())[1];
            int xmax = (box.high())[0];
            int ymax = (box.high())[1];
            int xsize = (xmax - xmin) + 1;
            int ysize = (ymax - ymin) + 1;

            FILE* fp = stderr;
            bool isFile = !file.empty();
            if (isFile) {
                fp = fopen(file.c_str(), "w");
            }

            if (!isFile) {
                fprintf(fp, "%s={", name.c_str());
            }
            for (unsigned i = 0; i < size; i++) {
                fprintf(fp, "%lg", data[i]);
                if (i < size - 1) {
                    fprintf(fp, ",");
                }
            }
            if (!isFile) {
                fprintf(fp, "}\n");
            } else {
                fclose(fp);
            }
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        Space newSpace(const BoxData<T, C, D, E>& boxData) {
            return newSpace("", boxData);
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        void addComponents(vector<Iter>& iters, vector<pdfg::Range>& ranges, vector<int>& offsets) {
            int csizes[3] = {E, D, C};
            char citers[3] = {'e', 'd', 'c'};
            const char* consts[3] = {"E", "D", "C"};

            for (unsigned c = 0; c < 3; c++) {
                if (csizes[c] > 1) {
                    Iter iter(citers[c]);
                    iters.emplace_back(iter);
                    Math upper = getUpper(csizes[c], string(consts[c]));
                    pdfg::Range range = (Int(0) <= iter < upper); //Int(csizes[c]));
                    ranges.emplace_back(range);
                    offsets.emplace_back(0);
                }
            }
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

            addComponents<T,C,D,E>(iters, ranges, offsets);

            Box box = boxData.box();
            Point low = box.low();
            Point high = box.high();

            // Add dimensions in reverse order
            for (int d = DIM-1; d >= 0; d--) {
                Iter iter(_iterNames[d]);
                iters.emplace_back(iter);
                Math upper = getUpper(high[d], "N", NUMCELLS);  // Int(high[d]);
                pdfg::Range range = (Int(low[d]) <= iter <= upper);
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
            Access access(space.name(), tuple, '(', offsets);

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
            // Add dimensions in reverse order
            for (int d = DIM-1; d >= 0; d--) {
                Iter iter(_iterNames[d]);
                Math upper = getUpper(high[d], "N", NUMCELLS);
                pdfg::Range range = (Int(low[d]) <= iter <= upper);
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
                cerr << "ERROR: No functor found for function '" << name << "'!\n";
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

            // Create iteration and output spaces.
            Space ispace = newSpace("I" + name, domain);
            Space ospace = newSpace(outName, outData);

            Macro macro = newMacro<T,C,D,E>(name, ispace, ospace, fxn);
            Comp comp(macro.name(), ispace, call(macro));
            return comp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outname, const Stencil<T>& stencil,
                     const BoxData<T,C,D,E>& src, T scale, const BoxData<T,C,D,E>& dest) {
            // Make a copy of m_coefs and multiply by a_scale
            vector<T> coefs = stencil.coefs();
            unsigned size = coefs.size();
            if (scale != 1) {
                for (int ii = 0; ii < size; ii++) {
                    coefs[ii] *= scale;
                }
            }

            // convert Points in m_offsets to linear offsets in a_src's data
//            int offsetFactor[DIM];
//            offsetFactor[0] = 1;
//            for (int ii = 1; ii < DIM; ii++) {
//                offsetFactor[ii] = offsetFactor[ii-1] * src.box().size(ii-1);
//            }
//
//            vector<int> offsets(size, 0);
//            vector<Point> offpoints = stencil.offsets();
//            for (int ii = 0; ii < size; ii++) {
//                for (int dir = 0; dir < DIM; dir++) {
//                    offsets[ii] += offpoints[ii][dir]*offsetFactor[dir];
//                }
//            }

            // apply the stencil
            Point srcRefRatio = stencil.srcRatio();
            Point destRefRatio = stencil.destRatio();
            Point destShift = stencil.destShift();
            vector<Point> offsets = stencil.offsets();

            // Get the source data space.
            auto itr = _dataPtrs.find(src.data());
            Space sspace;
            if (itr != _dataPtrs.end()) {
                sspace = getSpace(itr->second);
            } else {
                sspace = newSpace("", src);
            }

            // Get the destination space.
            Space dspace = newSpace(outname, dest);

            // Get the iteration space.
            Box box = stencil.indexRange(src.box());
            Space ispace = newSpace("I" + name, box);

            // Merge the iteration space w/ dest space as stencils need component loops.
            ispace.merge(dspace);
            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            // Calculate the number of dimensions in the component tensors.
            unsigned nCompDims = nComponentDims<T,C,D,E>();
            for (unsigned dim = nCompDims; dim < iters.size(); dim++) {
                int shift = destShift[dim - nCompDims];
                if (shift < 0) {
                    tuple[dim] = (tuple[dim] - (-shift));
                } else if (shift > 0) {
                    tuple[dim] = (tuple[dim] + shift);
                }
            }

            Access lhs(dspace.name(), tuple);

            Math rhs;
            for (int jj = 0; jj < size; jj++) {
                Real coef(coefs[jj]);
                Point point = offsets[jj] * srcRefRatio;
                //int offset = offsets[jj];

                tuple.clear();
                int n = 0;
                for (; n < nCompDims; n++) {
                    tuple.push_back(iters[n]);
                }
                for (int d = DIM - 1; n < iters.size() && d >= 0; n++, d--) {
                    int offset = point[d];
                    if (offset < 0) {
                        tuple.push_back(iters[n] - (-offset));
                    } else {
                        tuple.push_back(iters[n] + offset);
                    }
                }

                Access acc(sspace.name(), tuple);
                rhs = rhs + paren(coef * acc);
            }

            Math expr = (lhs = rhs);
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

            Access readAcc(rspace.name(), tuple);
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
                     const BoxData<T,C,D,E>& lhs, T rhs) {
            Box domain = lhs.box();
            Space lspace = getSpace(_dataPtrs[lhs.data()]);

            Space ispace = newSpace("I" + name, domain);
            ispace.merge(lspace);

            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            Access writeAcc(lspace.name(), tuple);
            Real scalar(rhs);
            Math expr(writeAcc, scalar, oper);

            Comp comp(name, ispace, (expr));
            return comp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const string& oper,
                     const BoxData<T,C,D,E>& lhs, const BoxData<T,C,D,E>& rhs) {
            Box lbox = lhs.box();
            Box rbox = rhs.box();
            Box domain;
            if (lbox == rbox) {
                domain = lbox;
            } else {
                domain = lbox & rbox;
            }

            Space lspace = getSpace(_dataPtrs[lhs.data()]);
            Space rspace = getSpace(_dataPtrs[rhs.data()]);

            return newComp<T,C,D,E>(name, outName, oper, domain, lspace, rspace);
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& stencilName, const string& stencilOut, const Stencil<T>& stencil, T scale,
                     const string& opName, const string& outName, const string& oper,
                     const BoxData<T,C,D,E>& lhs, const BoxData<T,C,D,E>& rhs) {
//            Box lbox = lhs.box();
//            Box rbox = rhs.box();
//            Box domain;
//            if (lbox == rbox) {
//                domain = lbox;
//            } else {
//                domain = lbox & rbox;
//            }
//            Space tspace = newSpace(stencilOut, lhs);

            // 1) Apply stencil computation with temporary intermediate output.
//            Comp newComp(const string& name, const string& outname, const Stencil<T>& stencil,
//                         const BoxData<T,C,D,E>& src, T scale, const BoxData<T,C,D,E>& dest)
            Comp stencilComp = newComp<T,C,D,E>(stencilName, stencilOut, stencil, rhs, scale, lhs);

            // 2) Apply operation...
            Comp operComp = newComp<T,C,D,E>(opName, outName, oper, lhs, rhs);

            return operComp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& outName, const string& oper,
                     Box domain, Space lspace, Space rspace) {
            Space ispace = newSpace("I" + name, domain);
            ispace.merge(lspace);

            vector<Iter> iters = ispace.iterators();
            vector<Expr> tuple(iters.begin(), iters.end());

            Access writeAcc(lspace.name(), tuple);
            Access readAcc(rspace.name(), tuple);
            Math expr(writeAcc, readAcc, oper);

            Comp comp(name, ispace, (expr));
            return comp;
        }

    protected:
        DataFlowFactory() = default;

        Math getUpper(int bound, const string& constName, int constValue = 0) {
            if (constValue == 0) {
                constValue = bound;
            }
            Math upper;
            Const con = getConstant(constName);
            int diff = bound - constValue;
            if (diff < 0) {
                upper = con - Int(-diff);
            } else {
                upper = con + Int(diff);
            }
            return upper;
        }

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

        template <typename T, unsigned int C, unsigned char D, unsigned char E>
        inline unsigned nComponentDims() {
            return ((E > 1) ? 1 : 0) + ((D > 1) ? 1 : 0) + ((C > 1) ? 1 : 0);
        }

#if DIM>2
        const char _iterNames[9] = {'x', 'y', 'z', 'i', 'j', 'k', 'm', 'n', 'p'};
#elif DIM>1
        const char _iterNames[6] = {'x', 'y', 'i', 'j', 'm', 'n'};
#else
        const char _iterNames[3] = {'x', 'i', 'm'};
#endif
        map<const void*, string> _dataPtrs;
        vector<string> _inputNames;
        size_t _inputPos;
    };
}

#endif //PROTO_PROTO_DATAFLOW_H
