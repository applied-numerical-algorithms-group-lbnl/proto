HDF5Handler::HDF5Handler()
{
    // define custom types for Proto::Point and Proto::Box
    char indices[] = "ijklmnopqrstuvwxyz";
    
    H5T_PROTO_POINT = H5Tcreate(H5T_COMPOUND, sizeof(Point));
    for (int ii = 0; ii < DIM; ii++)
    {
        std::string pvar = "intvect";
        pvar += indices[ii];
        H5Tinsert(H5T_PROTO_POINT, pvar.c_str(), HOFFSET(Point, m_tuple)+sizeof(int)*ii, H5T_NATIVE_INT); 
    }

    H5T_PROTO_BOX = H5Tcreate(H5T_COMPOUND, sizeof(Box));
    for (int ii = 0; ii < DIM; ii++)
    {
        std::string lowName = "lo_";
        lowName += indices[ii];
        H5Tinsert(H5T_PROTO_BOX, lowName.c_str(), HOFFSET(Box, m_low)+ii*sizeof(int), H5T_NATIVE_INT);
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        std::string highName = "hi_";
        highName += indices[ii];
        H5Tinsert(H5T_PROTO_BOX, highName.c_str(), HOFFSET(Box, m_high)+ii*sizeof(int), H5T_NATIVE_INT);
    }

    H5T_PROTO_STRING = H5Tcopy(H5T_C_S1);
    H5Tset_size(H5T_PROTO_STRING, 100);
    //H5T_PROTO_STRING = H5Tcreate(H5T_STRING, H5T_VARIABLE);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
void HDF5Handler::writeLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
        array<double, DIM>& a_dx,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    writeLevel(a_data, varnames, a_dx, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR, typename... Args>
void HDF5Handler::writeLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
        std::vector<std::string> a_varNames,
        array<double, DIM>& a_dx,
        std::string a_filename,
        Args... a_params)
{
    // PARSE THE FILE NAME
    char fname[100];
    if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
    {
        a_filename += ".hdf5";
    }
    
    // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma diagnostic ignored "-Wformat-security"
    sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop
   
    // CREATE THE FILE
#ifdef PR_MPI
    
    auto mpi_comm = MPI_COMM_WORLD;
    auto mpi_info = MPI_INFO_NULL;

    auto p_access = H5Pcreate(H5P_FILE_ACCESS);
    assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, p_access);
    assert(H5Pclose(p_access) >= 0);
#else
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif
    auto s_scalar = H5Screate(H5S_SCALAR);
    auto g_root = H5Gopen(file, "/", H5P_DEFAULT);
    
    // ATTRIBUTES
    auto att_filetype =       H5Acreate2(g_root, "filetype", H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_num_components = H5Acreate2(g_root, "num_components", H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_num_levels =     H5Acreate2(g_root, "num_levels", H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);

    std::string filetype = "VanillaAMRFileType";
    int numComponents = C;
    int numLevels = 1;
     
    assert(H5Awrite(att_filetype, H5T_PROTO_STRING, filetype.c_str()) >= 0);
    assert(H5Awrite(att_num_components, H5T_NATIVE_INT, &numComponents) >= 0);
    assert(H5Awrite(att_num_levels, H5T_NATIVE_INT, &numLevels) >= 0);

    assert(H5Aclose(att_filetype) >= 0);
    assert(H5Aclose(att_num_components) >= 0);
    assert(H5Aclose(att_num_levels) >= 0);
    
    while (a_varNames.size() < C)
    {
        int index = a_varNames.size() + 1;
        a_varNames.push_back("var" + to_string(index));
    }
    
    for (int ii = 0; ii < C; ii++)
    {
        char compName[100];
        sprintf(compName, "component_%i", ii);

        auto att_component = H5Acreate2(g_root, compName, H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
        std::string c = a_varNames[ii];
        assert(H5Awrite(att_component, H5T_PROTO_STRING, c.c_str()) >= 0);
        assert(H5Aclose(att_component) >= 0);
    }

    assert(H5Gclose(g_root) >= 0);
    
    // CHOMBO GROUP
    char g_chombo_name[] = "/Chombo_global";
    auto g_chombo = H5Gcreate(file, g_chombo_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto att_SpaceDim = H5Acreate2(g_chombo, "SpaceDim", H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_testReal = H5Acreate2(g_chombo, "testReal", H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);

    int SpaceDim = DIM;
    assert(H5Awrite(att_SpaceDim, H5T_NATIVE_INT, &SpaceDim) >= 0);

    assert(H5Aclose(att_SpaceDim) >= 0);
    assert(H5Aclose(att_testReal) >= 0);
    assert(H5Gclose(g_chombo) >= 0);
    // LEVEL_0
    addLevel(&file, a_data, a_dx, 0);
    
    assert(H5Sclose(s_scalar) >= 0);
    assert(H5Fclose(file) >= 0);
}
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void HDF5Handler::addLevel(hid_t* a_file,
        LevelBoxData<T, C, MEM, CTR>& a_data,
        array<double, DIM>& a_dx0,              // dx maybe ought to be type T instead of double
        int a_level)
{
    unsigned int numPatches = a_data.layout().size();
    unsigned int patchSize = a_data.patchSize();
    unsigned int numData = numPatches*patchSize;

    unsigned int numPatches_local = a_data.layout().mySize();
    unsigned int numData_local = numPatches_local*patchSize;

    auto s_scalar = H5Screate(H5S_SCALAR);
    
    // CREATE GROUP
    std::string g_level_name = "/level_" + to_string(a_level);
    auto g_level = H5Gcreate(*a_file, g_level_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // CREATE ATTRIBUTES
    auto att_dt   = H5Acreate2(g_level, "dt",   H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_time = H5Acreate2(g_level, "time", H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_dx   = H5Acreate2(g_level, "dx",   H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT); //FIXME: anisotropic
    auto att_ref_ratio   = H5Acreate2(g_level, "ref_ratio",   H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT); //FIXME: anisotropic
    auto att_prob_domain = H5Acreate2(g_level, "prob_domain", H5T_PROTO_BOX,  s_scalar, H5P_DEFAULT, H5P_DEFAULT); 

    double dt = 1.0;
    double time = 0.0;
    double dx = a_dx0[0];
    int refRatio = 1;
    Box problemDomain = a_data.layout().domain().box();

    assert(H5Awrite(att_dt,   H5T_NATIVE_DOUBLE, &dt) >= 0);
    assert(H5Awrite(att_time, H5T_NATIVE_DOUBLE, &time) >= 0);
    assert(H5Awrite(att_dx,   H5T_NATIVE_DOUBLE, &dx) >= 0);
    assert(H5Awrite(att_ref_ratio,   H5T_NATIVE_INT, &refRatio) >= 0);
    assert(H5Awrite(att_prob_domain, H5T_PROTO_BOX,  &problemDomain) >= 0);

    assert(H5Aclose(att_dt) >= 0);
    assert(H5Aclose(att_time) >= 0);
    assert(H5Aclose(att_dx) >= 0);
    assert(H5Aclose(att_ref_ratio) >= 0);
    assert(H5Aclose(att_prob_domain) >= 0);

    // CREATE DATASPACES
    hsize_t procDims[] = {numProc()};
    hsize_t offsetDims[] = {numPatches + 1};
    hsize_t boxDims[] = {numPatches};
    hsize_t dataDims[] = {numData};

    hsize_t boxDims_local[] = {numPatches_local};
    hsize_t dataDims_local[] = {numData_local};

    auto s_processors = H5Screate_simple(1, procDims, NULL);
    auto s_offsets    = H5Screate_simple(1, offsetDims, NULL);
    auto s_boxes      = H5Screate_simple(1, boxDims, NULL);
    auto s_data       = H5Screate_simple(1, dataDims, NULL);
    
    auto s_boxes_local = H5Screate_simple(1, boxDims_local, NULL);
    auto s_data_local  = H5Screate_simple(1, dataDims_local, NULL);

    // CREATE DATASETS
    auto ds_processors = H5Dcreate2(*a_file, (g_level_name + "/processors").c_str(),     H5T_NATIVE_INT, s_processors,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_boxes =      H5Dcreate2(*a_file, (g_level_name + "/boxes").c_str(),          H5T_PROTO_BOX, s_boxes,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_offsets =    H5Dcreate2(*a_file, (g_level_name + "/data:offsets=0").c_str(), H5T_NATIVE_LONG, s_offsets,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t ds_data;
    addDataSet<T>(&ds_data, a_file, (g_level_name + "/data:datatype=0").c_str(), &s_data);

    // GET RAW DATA
    long int offsetData[numPatches + 1];
    Box boxData[numPatches_local];
    T rawData[numData_local];
    
    offsetData[0] = 0;
    int offset = 0;
    if (true) {
        for (int ii = 1; ii <= numPatches; ii++)
        {
            offset += patchSize;
            offsetData[ii] = offset;
        }
        int ii = 0;
        for (auto iter = a_data.begin(); *iter != iter.end(); ++iter, ii++)
        {
            boxData[ii] = a_data[*iter].box().grow(-a_data.ghost()); //this should just be a_data.layout()[*iter]...
        }
    } else {
        int ii = 0;
        for (auto iter = a_data.begin(); *iter != iter.end(); ++iter, ii++)
        {
            offset += a_data.size(*iter)/sizeof(T);
            offsetData[ii + 1] = offset;

            boxData[ii] = a_data[*iter].box().grow(-a_data.ghost());
        }
    }
    
    a_data.linearOut(rawData);

    // WRITE DATA
    // FIXME: write processor data
    assert(H5Dwrite(ds_offsets, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, offsetData) >= 0);
    #ifdef PR_MPI
    int boxStart[] = {a_data.layout().offset()};
    int dataStart[] = {a_data.offset()};
    int stride[] = {1};
    
    auto slab_boxes = H5Dget_space(ds_boxes); 
    assert(H5Sselect_hyperslab(slab_boxes, H5S_SELECT_SET, boxStart,  stride, boxDims_local, NULL);
    auto slab_data = H5Dget_space(ds_data); 
    assert(H5Sselect_hyperslab(slab_data,  H5S_SELECT_SET, dataStart, stride, dataDims_local, NULL);

    assert(H5Dwrite(ds_boxes,   H5T_PROTO_BOX, s_boxes_local, slab_boxes, H5P_DEFAULT, boxData) >= 0);
    writeData(a_file, &ds_data, &s_data_local, &slab_data, rawData);

    assert(H5Sclose(slab_boxes) >= 0);
    assert(H5Sclose(slab_data) >= 0);
    #else
    assert(H5Dwrite(ds_boxes,   H5T_PROTO_BOX, H5S_ALL, H5S_ALL, H5P_DEFAULT, boxData) >= 0);
    writeData(a_file, &ds_data, &s_data, &s_data, rawData);
    #endif
    

    // CLEANUP
    assert(H5Dclose(ds_boxes) >= 0);
    assert(H5Dclose(ds_offsets) >= 0);
    assert(H5Dclose(ds_data) >= 0);

    assert(H5Sclose(s_offsets) >= 0);
    assert(H5Sclose(s_processors) >= 0);
    assert(H5Sclose(s_boxes) >= 0);
    assert(H5Sclose(s_data) >= 0);
    assert(H5Sclose(s_boxes_local) >= 0);
    assert(H5Sclose(s_data_local) >= 0);

    // METADATA
    std::string g_level_meta_name = g_level_name + "/data_attributes";
    auto g_level_meta = H5Gcreate(*a_file, g_level_meta_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto att_comps = H5Acreate2(g_level_meta, "comps", H5T_NATIVE_INT,  s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_ghost = H5Acreate2(g_level_meta, "ghost", H5T_PROTO_POINT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_object_type = H5Acreate2(g_level_meta, "object_type", H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_outputGhost = H5Acreate2(g_level_meta, "outputGhost", H5T_PROTO_POINT,  s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    
    int comps = C;
    Point ghost = a_data.ghost();
    Point outputGhost = a_data.ghost(); //FIXME No idea what this attribute is meant to be

    assert(H5Awrite(att_comps, H5T_NATIVE_INT, &comps) >= 0);
    assert(H5Awrite(att_ghost, H5T_PROTO_POINT, &ghost) >= 0);
    assert(H5Awrite(att_object_type, H5T_PROTO_STRING, "FArrayBox") >= 0);
    assert(H5Awrite(att_outputGhost, H5T_PROTO_POINT, &outputGhost) >= 0);

    assert(H5Aclose(att_comps) >= 0);
    assert(H5Aclose(att_ghost) >= 0);
    assert(H5Aclose(att_object_type) >= 0);
    assert(H5Aclose(att_outputGhost) >= 0);

    // FINAL CLEAN UP
    assert(H5Gclose(g_level_meta) >= 0);
    assert(H5Gclose(g_level) >= 0);
    assert(H5Sclose(s_scalar) >= 0);
}

  
  
    
/*
void HDF5Handler::writeData( hid_t* a_file,
        std::string a_dataName,
        hid_t* a_space,
        double* a_data)
{
    auto ds_data = H5Dcreate2(*a_file, a_dataName.c_str(), H5T_NATIVE_DOUBLE, *a_space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(H5Dwrite(ds_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data) >= 0);
    assert(H5Dclose(ds_data) >= 0);
}

void HDF5Handler::writeData( hid_t* a_file,
        std::string a_dataName,
        hid_t* a_space,
        int* a_data)
{
    auto ds_data = H5Dcreate2(*a_file, a_dataName.c_str(), H5T_NATIVE_INT, *a_space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(H5Dwrite(ds_data, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data) >= 0);
    assert(H5Dclose(ds_data) >= 0);
}
*/

template<typename T>
void
HDF5Handler::addDataSet(hid_t* a_dataSet,
        hid_t* a_file,
        std::string a_name,
        hid_t* a_dataSpace)
{
    PROTO_ASSERT(false, "HDF5Handler::addDataSet | Error: Unrecognized data type."); 
}

template<>
void
HDF5Handler::addDataSet<double>(hid_t* a_dataSet,
        hid_t* a_file,
        std::string a_name,
        hid_t* a_dataSpace)
{
    *a_dataSet = H5Dcreate2(*a_file, a_name.c_str(), H5T_NATIVE_DOUBLE, *a_dataSpace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

template<>
void
HDF5Handler::addDataSet<int>(hid_t* a_dataSet,
        hid_t* a_file,
        std::string a_name,
        hid_t* a_dataSpace)
{
    *a_dataSet = H5Dcreate2(*a_file, a_name.c_str(), H5T_NATIVE_INT, *a_dataSpace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

template<typename T>
void
HDF5Handler::writeData(hid_t* a_file,
        hid_t* a_dataSet,
        hid_t* a_domainSpace,
        hid_t* a_localSpace,
        T* a_data)
{
    PROTO_ASSERT(false, "HDF5Handler::writeData | Error: Unrecognized data type."); 
}

template<>
void
HDF5Handler::writeData<double>(hid_t* a_file,
        hid_t* a_dataSet,
        hid_t* a_domainSpace,
        hid_t* a_localSpace,
        double* a_data)
{
    assert(H5Dwrite(*a_dataSet, H5T_NATIVE_DOUBLE, *a_localSpace, *a_domainSpace, H5P_DEFAULT, a_data) >= 0);
}

template<>
void
HDF5Handler::writeData<int>(hid_t* a_file,
        hid_t* a_dataSet,
        hid_t* a_domainSpace,
        hid_t* a_localSpace,
        int* a_data)
{
    assert(H5Dwrite(*a_dataSet, H5T_NATIVE_INT, *a_localSpace, *a_domainSpace, H5P_DEFAULT, a_data) >= 0);
}
