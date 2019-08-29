
/// IMEX Integrator Data Holder
/**
    This pure abstract class documents the interface needed by a class to be a valid
    data holder for IMEXIntegrator. 

    The data holder NEED NOT INHERIT FROM THIS CLASS. This file is for documentation
    purposes only.
*/
class IMEXData {
    public:
    /// Memory Allocation
    /**
        Allocate space for this according to an existing instance. No data is copied.
        TODO: This should probably be renamed "allocate" since it doesn't copy data
    */    
    void define(IMEXData& a_data) = 0;

    /// Copy
    /**
        Copy the data in this to another IMEXData instance.
    */
    void copyTo(IMEXData& a_data) = 0;

    /// Increment
    /**
        Add another IMEXData instance to this with an optional scaling.
    */
    void increment(IMEXData& a_data, double a_scale = 1.0) = 0;
};

