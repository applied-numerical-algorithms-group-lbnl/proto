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

using std::string;
using std::ostringstream;

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
        Space create(const BoxData<T, C, D, E>& boxData, const string& name = "") {
            string sname = name;
            if (sname.empty()) {
                ostringstream os;
                os << boxData.data();
                sname = os.str();
            }

            Box box = boxData.box();
            Point low = box.low();
            Point high = box.high();

            vector<Constr> constrs;
            for (int d = 0; d < DIM; d++) {
                Iter iter(_iterNames[d]);
                pdfg::Range range = (Int(low[d]) <= iter <= Int(high[d]));
                constrs.emplace_back(range.lower());
                constrs.emplace_back(range.upper());
            }

            Space space(sname, constrs);
            return space;
        }

    protected:
        DataFlowFactory() = default;

        const char _iterNames[6] = {'i', 'j', 'k', 'l', 'm', 'n'};
    };
}

#endif //PROTO_PROTO_DATAFLOW_H
