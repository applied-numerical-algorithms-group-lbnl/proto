//
// Created by edavis on 5/31/19.
//

#pragma once
#ifndef PROTO_PROTO_DATAFLOW_H
#define PROTO_PROTO_DATAFLOW_H

#include <pdfg/GraphIL.hpp>

//using pdfg::Iter;
//using pdfg::Constr;
//using pdfg::Space;
//using pdfg::Comp;
using namespace pdfg;

using std::cerr;
using std::endl;
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
            Box box = boxData.box();
            return newSpace(boxData.box());
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
            return space;
        }

//        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1, typename Func, typename... Srcs>
//        Comp newComp(const Func& func, Srcs&... srcs) {
//            ostringstream os;
//            os << &func;
//            return newComp(os.str(), func, srcs...);
//        }

        template<typename T, unsigned int C, unsigned char D = 1, unsigned char E = 1, typename Func, typename... Srcs>
        Comp newComp(const string& name, const Func& func, Srcs&... srcs) {
            // Compute intersection to get iteration space for 'forall'.
            _scalars.clear();
            Box intersect;
            parseArgs(intersect, srcs...);
            Space ispace = newSpace("I" + name, intersect);

            // TODO: Make this much better!
            Pointer fptr((void*) func);
            Comp comp(name, ispace, (call(name, {fptr})));

            return comp;
        }

    protected:
        DataFlowFactory() = default;

        template <typename T, unsigned int C, unsigned char D, unsigned char E>
        void parse(Box& domain, const BoxData<T,C,D,E>& data) {
            // Recursively compute domain intersection...
            if (domain.empty()) {
                domain = data.box();
            } else {
                domain = domain & data.box();
            }
        }

        template <typename T>
        void parse(Box& a_domain, T scalar) {
            // Satisfy compiler, we can convert this into a 'Real' object, but then what?
            Real real(scalar);
            _scalars.push_back(real);
        }

        template <typename Last>
        void parseArgs(Box& domain, Last& last) {
            parse(domain, last);
        }

        template <typename T, unsigned int C, unsigned char D, unsigned char E, typename... Rest>
        void parseArgs(Box& domain, const BoxData<T,C,D,E>& first, Rest&... rest) {
            parse(domain, first);
            parseArgs(domain, rest...);
        }

        template <typename T, unsigned int C, unsigned char D, unsigned char E>
        void parseArgs(Box& domain, const BoxData<T,C,D,E>& first) {
            parse(domain, first);
        }

        template <typename First , typename... Rest>
        void parseArgs(Box& domain, First& first, Rest&... rest)
        {
            parse(domain, first);
            parseArgs(domain, rest...);
        }

        const char _iterNames[6] = {'i', 'j', 'k', 'l', 'm', 'n'};
        vector<Real> _scalars;
    };
}

#endif //PROTO_PROTO_DATAFLOW_H
