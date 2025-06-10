
template<typename MAP, MemType MEM>
template <typename T, unsigned int C, MemType DATAMEM, Centering CTR>
void MBLevelArgs<MAP,MEM>::Set(const std::string& label, const std::shared_ptr<MBLevelBoxData<T, C, DATAMEM, CTR>>& inData) {
    PROTO_ASSERT(layout().compatible(inData->layout()),
        "MBLevelArgs::Set | Error: Level data is not compatible with map layout");
    _data[label] = std::make_shared<MBLevelBoxDataWrapper<T, C, DATAMEM, CTR>>(inData);
}
template<typename MAP, MemType MEM>
template<typename T, unsigned int C, MemType DATAMEM, Centering CTR>
MBLevelBoxData<T,C,DATAMEM,CTR>& MBLevelArgs<MAP,MEM>::Get(const std::string& label) {
    auto iter = _data.find(label);
    PROTO_ASSERT(iter != _data.end(),
        "MBLevelArgs::Get | Error: Could not find data label: %s", label.c_str());
    auto wrapper = std::dynamic_pointer_cast<MBLevelBoxDataWrapper<T,C,DATAMEM,CTR>>(iter->second);
    PROTO_ASSERT(wrapper,
        "MBLevelArgs::Get | Error: Data type mismatch for label: %s", label.c_str());
    return *(wrapper->data);
}
template<typename MAP, MemType MEM>
template <typename T>
void MBLevelArgs<MAP,MEM>::Set(const std::string& label, const T& value) {
    constants[label] = value;
}
template<typename MAP, MemType MEM>
template<typename T>
T MBLevelArgs<MAP,MEM>::Get(const std::string& label) const {
    auto iter = constants.find(label);
    PROTO_ASSERT(iter != constants.end(),
        "MBLevelArgs::Get | Error: Could not find constant label: %s", label.c_str());
    return std::any_cast<T>(iter->second);
}
