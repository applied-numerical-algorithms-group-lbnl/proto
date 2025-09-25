#include "MBPoisson.H"
#include "Proto.H"
#include "TestFunctions.H"
#include "BoxOp_MBLaplace.H"
//#include "MBMap_Identity.H"



template <typename T, unsigned int C, MemType MEM, Centering CTR, typename MAP>
void initializeData(
        MBLevelBoxData<T, C, MEM, CTR> &phi,
        MBLevelBoxData<T, C, MEM, CTR> &lphi,
        MBLevelMap<MAP, MEM>& map, Point ghost)
{
        auto& layout = phi.layout();
        double dx = 0;
        Array<double, DIM> offset{dx, dx ,0,0,0,0};
        auto C2C = Stencil<double>::CornersToCells(4);
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto &phi_i = phi[iter];
            auto &rhs_i = lphi[iter];
            Box b_i = C2C.domain(layout[iter]).grow(ghost);
            BoxData<double, DIM> x_i(b_i.grow(PR_NODE));
            BoxData<double, 1> J_i(x_i.box());
            // compute coordinate / jacobian values from map
            map.apply(x_i, J_i, block);
            BoxData<double, 1> phi_node = forall<double, 1>(f_bell, x_i, offset);
            BoxData<double, 1> lphi_node = forall<double, 1>(f_Lbell, x_i, offset);

            phi_i |= C2C(phi_node);
            rhs_i |= C2C(lphi_node);
        }
}

int main(int argc, char* argv[])
{
    #if PR_VERBOSE > 0
    HDF5Handler h5;
    #endif

    int domainSize = 16;
    int boxSize = 8;
    constexpr int numBlocks = 4;
    //constexpr int numBlocks = 5;
    int numLevels = log(domainSize)/log(2.0)-2;
    Point refRatio = Point::Ones(2);
   typedef BoxOp_MBLaplace<double, MBMap_XPointRigid<numBlocks, HOST>> OP;
    auto domain = XPoint<numBlocks>::Domain(domainSize, boxSize);
    MBDisjointBoxLayout layout(domain, Point::Ones(boxSize));

    MBMultigrid<BoxOp_MBLaplace, MBMap_XPointRigid<numBlocks, HOST>, double> mg(layout, refRatio, numLevels);
    MBLevelBoxData<double, 1, HOST> phi(layout, OP::ghost());
    MBLevelBoxData<double, 1, HOST> rhs(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> res(layout, Point::Zeros());

    MBLevelMap<MBMap_XPointRigid<numBlocks, HOST>, HOST> map(layout, OP::ghost());

    initializeData(phi, rhs, map, OP::ghost());
    
    mg.op(numLevels-1)(rhs, phi);
    phi.setVal(0);
    res.setVal(0);
#if PR_VERBOSE > 0
    h5.writeMBLevel({"rhs"}, map, rhs, "RHS_LAPLACE_XPOINT");
    h5.writeMBLevel({"phi"}, map, phi, "PHI_LAPLACE_XPOINT");
    h5.writeMBLevel({"residual"}, map, res, "RES_LAPLACE_XPOINT");
#endif
    mg.solve(phi, rhs, 20, 1e-10);
}
