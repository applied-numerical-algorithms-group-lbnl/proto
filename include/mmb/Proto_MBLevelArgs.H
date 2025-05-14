#pragma once
#ifndef _PROTO_MB_LEVEL_ARGS_
#define _PROTO_MB_LEVEL_ARGS_

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

    template<typename MAP, MemType MEM>
    class MBLevelArgs {
    private:

        std::shared_ptr<MBLevelMap<MAP,MEM>> _map;
        std::unordered_map<std::string, std::shared_ptr<DataWrapper>> _data;
        std::unordered_map<std::string, std::any> constants;
    public:
        
        MBLevelArgs(std::shared_ptr<MBLevelMap<MAP,MEM>>& map) { _map = map; }

        MBLevelArgs(const MBLevelArgs& rhs) = default;

        inline const MBDisjointBoxLayout& layout() { return _map->layout(); }
        inline const MBLevelMap<MAP,MEM>& map() { return *_map; }
        template <typename T, unsigned int C, MemType DATAMEM, Centering CTR=PR_CELL>
        inline void Set(const std::string& label, const std::shared_ptr<MBLevelBoxData<T, C, DATAMEM, CTR>>& inData);

        template<typename T, unsigned int C, MemType DATAMEM, Centering CTR=PR_CELL>
        inline MBLevelBoxData<T,C,DATAMEM,CTR>& Get(const std::string& label); 

        template <typename T>
        inline void Set(const std::string& label, const T& value); 

        template<typename T>
        inline T Get(const std::string& label) const;
    };

    #include "implem/Proto_MBLevelArgsImplem.H"
}
#endif