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
//#include <initializer_list>
using std::initializer_list;

namespace Proto {
    class DataFlowFactory {
    public:
        static DataFlowFactory& get() {
            static DataFlowFactory instance; // Guaranteed to be destroyed.
            return instance;            // Instantiated on first use.
        }

        DataFlowFactory(DataFlowFactory const&) = delete;         // Don't Implement

        void operator=(DataFlowFactory const&) = delete;     // Don't implement

        void init(const string& name) {
            pdfg::init(name);
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
            return newSpace(sname, boxData.box());
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

            // Compute intersection to get iteration space for 'forall'.
            Box intersect;
            pdfg::Func fxn(name);
            parseArgs(intersect, fxn, srcs...);

            if (domain.empty()) {
                domain = intersect;
            }

            BoxData<T,C,D,E> result(domain);
            Space ispace = newSpace("I" + name, domain);
            Space ospace = newSpace(outName, outData);
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

            Space A = dspace;
            Iter i = iters[0];
            Iter j = iters[1];
            Math expect = (A(i,j) = (A(i,j) + A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i,j-1)) * 0.2);

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
            if (scale == 1.0) {
                scale /= (T) coefs.size();
            }

            Math expr = lhs = rhs * scale;
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

        const char _iterNames[6] = {'i', 'j', 'k', 'l', 'm', 'n'};
        map<const void*, string> _dataPtrs;
        vector<string> _inputNames;
        size_t _inputPos;
    };
}

#endif //PROTO_PROTO_DATAFLOW_H
