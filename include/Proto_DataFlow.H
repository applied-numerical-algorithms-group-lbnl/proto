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
            char str[128];
            sprintf(str, "%p", boxData.data());
            return newSpace(string(str), boxData);
        }

        template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
        Space newSpace(const string& name, const BoxData<T, C, D, E>& boxData) {
            _dataPtrs[boxData.data()] = name;
            Box box = boxData.box();
            return newSpace(name, boxData.box());
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
        Comp newComp(const string& name, const string& output, const Func& func, Srcs&... srcs) {
            // Compute intersection to get iteration space for 'forall'.
            Box domain;
            pdfg::Func fxn(name);

            parseArgs(domain, fxn, srcs...);
            BoxData<T,C,D,E> result(domain);
            Space ispace = newSpace("I" + name, domain);
            Space ospace = newSpace(output, domain);
            fxn.eval();

            Comp comp(name, ispace, call(ospace, fxn));
            return comp;
        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1>
        Comp newComp(const string& name, const string& output, const Stencil<T>& stencil, const BoxData<T, C, D, E>& data) {
            auto itr = _dataPtrs.find(data.data());
            Space dspace;
            if (itr != _dataPtrs.end()) {
                dspace = getSpace(itr->second);
            } else {
                dspace = newSpace(data);
            }

            Box domain = data.box();
            Space ispace = newSpace("I" + name, domain);
            vector<Iter> iters = ispace.iterators();

            Space A = dspace;
            Iter i = iters[0];
            Iter j = iters[1];
            Math expect = (A(i,j) = (A(i,j) + A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i,j-1)) * 0.2);

            vector<Expr> tuple(iters.begin(), iters.end());

            // TODO: Calculate output space size for stencil...
            Space ospace = newSpace(output, domain);
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
            Real scale(1.0 / (T) coefs.size());
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
            // Add data space to function as arg...
            string name = _dataPtrs[data.data()];
            Space space = getSpace(name);
            func.add(space);
            int stop = 1;
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
    };
}

#endif //PROTO_PROTO_DATAFLOW_H
