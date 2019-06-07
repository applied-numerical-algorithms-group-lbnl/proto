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

// Macros (for later)
#define consToPrim(W,U,gamma,y,x) {\
(W)(0,(y),(x))=(U)(0,(y),(x));\
(W)(1,(y),(x))=(U)(1,(y),(x)) / (U)(0,(y),(x));\
(W)(2,(y),(x))=(U)(2,(y),(x)) / (U)(0,(y),(x));\
(W)(3,(y),(x))=((U)(3(y),(x))-0.5*(U)(0,(y),(x))*(((W)(1,(y),(x))*(W)(1,(y),(x)))+((W)(2,(y),(x))*(W)(2,(y),(x)))))*((gamma)-1.0);\
}

namespace Proto {
    class DataFlowFactory {
    public:
        static DataFlowFactory& get() {
            static DataFlowFactory instance; // Guaranteed to be destroyed.
            return instance;            // Instantiated on first use.
        }

        DataFlowFactory(DataFlowFactory const&) = delete;         // Don't Implement

        void operator=(DataFlowFactory const&) = delete;     // Don't implement

        void init(const string& name, const string& retname = "", const string& datatype = "", const string& indextype = "") {
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

#define offset2(i,j,M) ((j)+(i)*(M))
#define offset3(i,j,k,M,N) ((k)+((j)+(i)*(M))*(N))
#define offset4(i,j,k,l,M,N,P) ((l)+((k)+((j)+(i)*(M))*(N))*(P))

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

            // Create output space.
            Space ospace = newSpace(outName, outData);

            // Compute intersection to get iteration space for 'forall'.
            Box intersect;
            pdfg::Func fxn(name);
            //fxn.add(ospace);
            parseArgs(intersect, fxn, srcs...);

            if (domain.empty()) {
                domain = intersect;
            }

            BoxData<T,C,D,E> result(domain);
            Space ispace = newSpace("I" + name, domain);
            for (const auto& iter : ispace.iterators()) {
                fxn.add(iter);
            }
            fxn.eval();

            Comp comp(name, ispace, call(ospace, fxn));
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

            Access writeAcc(lspace, {Int(0)});
            Access readAcc(rspace, tuple);
            Math expr(writeAcc, readAcc, oper);

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
