#pragma once
#ifndef _PROTO_LEVEL_ARGS_
#define _PROTO_LEVEL_ARGS_

#include <string>
#include <any>
#include <unordered_map>

#include "Proto.H"

namespace Proto
{
    struct DataWrapper {
        virtual ~DataWrapper() = default;
    };

    template <typename T, unsigned int C, MemType MEM, Centering CTR=PR_CELL>
    struct MBLevelBoxDataWrapper : public DataWrapper {
        
        std::shared_ptr<MBLevelBoxData<T,C,MEM,CTR>> data;
        
        explicit MBLevelBoxDataWrapper(std::shared_ptr<MBLevelBoxData<T,C,MEM,CTR>> d)
        : data(std::move(d)) {}
    };

    template<template<MemType> typename MAP, MemType MEM>
    class MBLevelArgs {
    private:

        std::shared_ptr<MBLevelMap<MAP,MEM>> _map;
        std::unordered_map<std::string, std::shared_ptr<DataWrapper>> _data;
        std::unordered_map<std::string, std::any> constants;
    public:
        
        MBLevelArgs(std::shared_ptr<MBLevelMap<MAP,MEM>>& map)
        {
            _map = map;
        }

        const MBDisjointBoxLayout& layout() { return _map->layout(); }
        const MBLevelMap<MAP,MEM>& map() { return *_map; }
        template <typename T, unsigned int C, MemType DATAMEM, Centering CTR=PR_CELL>
        void Set(const std::string& label, const std::shared_ptr<MBLevelBoxData<T, C, DATAMEM, CTR>>& inData) {
            PROTO_ASSERT(layout().compatible(inData->layout()),
                "MBLevelArgs::Set | Error: Level data is not compatible with map layout");
            _data[label] = std::make_shared<MBLevelBoxDataWrapper<T, C, DATAMEM, CTR>>(inData);
        }

        template<typename T, unsigned int C, MemType DATAMEM, Centering CTR=PR_CELL>
        MBLevelBoxData<T,C,DATAMEM,CTR>& Get(const std::string& label) const {
            auto iter = _data.find(label);
            PROTO_ASSERT(iter != _data.end(),
                "MBLevelArgs::Get | Error: Could not find data label: %s", label.c_str());
            auto wrapper = std::dynamic_pointer_cast<MBLevelBoxDataWrapper<T,C,DATAMEM,CTR>>(iter->second);
            PROTO_ASSERT(wrapper,
                "MBLevelArgs::Get | Error: Data type mismatch for label: %s", label.c_str());
            return *(wrapper->data);
        }

        template <typename T>
        void Set(const std::string& label, const T& value) {
            constants[label] = value;
        }

        template<typename T>
        T Get(const std::string& label) const {
            auto iter = constants.find(label);
            PROTO_ASSERT(iter != constants.end(),
                "MBLevelArgs::Get | Error: Could not find constant label: %s", label.c_str());
            return std::any_cast<T>(iter->second);
        }
    };

}
#endif