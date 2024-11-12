#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

namespace {

    
    enum BCType {
        DIRICHLET,
        NEUMANN
    };

    // Data holder on a single face
    template<typename T, unsigned int C, MemType MEM>
    struct BC
    {
        BC(Face a_face, BCType a_type, Box a_interior) : face(a_face), type(a_type)
        {
            Point dir = Point::Basis(face.dir, face.sign());
            Box b = 
                face.sign() < 0 ? a_interior.edge(dir, 1)
              : a_interior.adjacent(dir, 1);
            data.define(b);
        }

        const Face face;        // which boundary is this
        const BCType type;      // what kind of data is this
        BoxData<T,C,MEM> data;  // the data
    };
    
    // Data holders for a patch (multiple faces / BC types possible)
    template<typename T, unsigned int C, MemType MEM>
    struct BCSet
    {
        public:

            BCSet(Box a_interior) : m_interior(a_interior) {}

            void add(Face a_face, BCType a_type)
            {
                m_data[faceToIndex(a_face)] = std::make_shared<BC<T,C,MEM>>(a_face, a_type, m_interior);
            }
            bool contains(Face a_face)
            {
                return m_data[faceToIndex(a_face)] != nullptr;
            }
            BC<T,C,MEM>& get(Face a_face)
            {
                return *m_data[faceToIndex(a_face)];
            }
            Box interior() const { return m_interior; }
            void print()
            {
                for (int dd = 0; dd < DIM; dd++)
                {
                    Face L(dd,Side::Lo);
                    pr_out() << "Dir: " << dd << " | Low Side" << std::endl;
                    if (contains(L))
                    {
                        get(L).data.printData();
                    } else {
                        pr_out() << "NO DATA FOUND" << std::endl;
                    }
                    Face H(dd,Side::Hi);
                    pr_out() << "Dir: " << dd << " | High Side" << std::endl;
                    if (contains(H))
                    {
                        get(H).data.printData();
                    } else {
                        pr_out() << "NO DATA FOUND" << std::endl;
                    }
                }
            }
        private:
            int faceToIndex(Face a_face)
            {
                return a_face.sign() < 0 ? a_face.dir : a_face.dir + DIM;
            }
        
            const Box m_interior;
            Array<std::shared_ptr<BC<T,C,MEM>>, 2*DIM> m_data;
    };

   
#if DIM == 3
    template<typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
    BoxData<T,C,MEM,D,E> Product(
            const BoxData<T,DIM,MEM>& a_L4,
            const BoxData<T,DIM,MEM>& a_R4,
            const BoxData<T,DIM,MEM>& a_L2,
            const BoxData<T,DIM,MEM>& a_R2,
            int                       a_edgeDir)
    {
    }

    template<typename T, MemType MEM>
    BoxData<T,DIM,MEM> CrossProduct(
            const BoxData<T,DIM,MEM>& a_L4,
            const BoxData<T,DIM,MEM>& a_R4,
            const BoxData<T,DIM,MEM>& a_L2,
            const BoxData<T,DIM,MEM>& a_R2,
            int                       a_edgeDir)
    {
    }
#endif


    template<typename T, MemType MEM>
    void Cofactor(BoxData<T,DIM,MEM>& a_NT, const BoxData<T,DIM,MEM>& a_X, int a_dir)
    {
#if DIM==2
        int d_tan = (a_dir + 1) % DIM;
        auto D = Stencil<T>::FluxDivergence(d_tan);
        for (int di = 0; di < DIM; di++)
        {
            int dj = (di + 1) % DIM;
            auto NTi = slice(a_NT,di);
            auto Xj = slice(a_X,dj);
            T sign = (dj == d_tan) ? 1.0 : -1.0;
            //NTi = dXj/dXn
            if (di == 0)
            {
                NTi |= D(Xj, sign);
            } else {
                NTi += D(Xj, sign);
            }
        }
#elif DIM==3
#endif
    }
    
    template<typename T, MemType MEM>
    BoxData<T,DIM,MEM> Cofactor(const BoxData<T,DIM,MEM>& a_X, int a_dir)
    {
#if DIM==2
        int d_tan = (a_dir + 1) % DIM;
        auto D = Stencil<T>::FluxDivergence(d_tan);
        BoxData<T,DIM,MEM> NT(D.range(a_X.box()));
        Cofactor(NT,a_X,a_dir);
        return NT;
#elif DIM==3
#endif

    }
    
    template<typename T, MemType MEM>
    FluxBoxData<T,DIM,MEM> Cofactor(const BoxData<T,DIM,MEM>& a_X)
    {
#if DIM==2
        Box B = a_X.box().extrude(Point::Ones(),-2);
        FluxBoxData<T,DIM,MEM> NT(B);
        for (int dd = 0; dd < DIM; dd++)
        {
            Cofactor(NT[dd], dd);
        }
        return NT;
#elif DIM==3
#endif
    }

    template<typename T, MemType MEM>
    BoxData<T,DIM,MEM,DIM> CofactorMatrix(
            const FluxBoxData<T,DIM,MEM>& a_NTFace,
            int a_norm)
    {
        Array<Stencil<T>,DIM> S;
        
        for (int dd = 0; dd < DIM; dd++)
        {
            int di = (a_norm + dd) % DIM;
            if (di == a_norm)
            {
                S[di] = Stencil<T>::Identity();
            } else {
                S[di] = Stencil<T>::CellToFace(a_norm) * Stencil<T>::faceToCell(di, 4);
            }
            std::cout << "n: " << a_norm << " | d: " << dd << " | span: " << S[di].span() << std::endl;
        }
        
        Box B = a_NTFace.box().grow(1);
        for (int dd = 0; dd < DIM; dd++)
        {
            B &= S[dd].range(a_NTFace[dd].box());
        }

        BoxData<T,DIM,MEM,DIM> NTOut(B);
        
        for (int dd = 0; dd < DIM; dd++)
        {
            auto NTOut_d = plane(NTOut,dd);
            NTOut_d |= S[dd](a_NTFace[dd]);
        }
        
        return NTOut;
    }

    template<typename T, unsigned int C, MemType MEM>
    void Gradient(
            FluxBoxData<T,C,MEM>&               a_grad,
            const BoxData<T,C,MEM>&             a_state,
            const BCSet<T,C,MEM>&               a_bcs,
            const BoxData<double, DIM, MEM>&    a_X,
            const BoxData<double, 1, MEM>&      a_J,
            Array<T,DIM>                        a_dx)
    {
        //===========================================================================
        //  Compute Metrics

        FluxBoxData<T,DIM,MEM> NT(a_grad.box());
        //FluxBoxData<T,1,MEM> Jf(a_grad.box());
        //FluxBoxData<T,DIM,MEM,DIM> NTMatrix(a_grad.box());
        std::vector<BoxData<T,DIM,MEM,DIM>> NTMatrix;
        std::vector<BoxData<T,1,MEM>> Jf;

        for (int dir = 0; dir < DIM; dir++)
        {
            auto C2F = Stencil<double>::CellToFace(dir);
            auto NT_dir = Operator::cofactor(a_X, dir, a_dx); 
            NT_dir.copyTo(NT[dir]);
            Jf.push_back(C2F(a_J));
        }
        for (int dir = 0; dir < DIM; dir++)
        {
            NTMatrix.push_back(Operator::cofactorMatrix(NT, dir));
        }

        //===========================================================================
        //  Compute Interior Gradient

        FluxBoxData<T,DIM,MEM> flux(a_grad.box());
#if 1
        for (int dir = 0; dir < DIM; dir++)
        {
            auto faceGrad =  Operator::_faceGradxPhi(
                    a_state, a_state, NTMatrix[dir], NTMatrix[dir],
                    Jf[dir], Jf[dir], dir);
            auto tmp = Operator::_faceMatrixProductAB(
                    faceGrad,NT[dir],faceGrad,NT[dir],dir);
            tmp.copyTo(a_grad[dir]);
            a_grad[dir] /= a_dx[dir];
        }
#else
#endif

        //===========================================================================
        //  Apply Boundary Conditions

        //===========================================================================
        //  Compute Product with NT
        
    }

}

TEST(BCDemo, Demo) {
    constexpr int C = 2;
    int domainSize = 8;
    int boxSize = 8;
    int numGhost = 2;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    
    Point ghost = Point::Ones(numGhost);
    auto domain = buildXPoint(domainSize);
    Point boxSizeV = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeV);

    MBLevelMap<MBMap_XPointRigid, HOST> map;
    map.define(layout, ghost);

    MBLevelBoxData<double, C, HOST> phi(layout, ghost); 
    MBLevelBoxData<double, C, HOST> lphi(layout, Point::Zeros()); 
    
    //------------------------------------------------------------------------------
    // Initialize Data
    auto C2C = Stencil<double>::CornersToCells(4);
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& phi_i = phi[iter];
        Box b_i = C2C.domain(layout[iter]).grow(numGhost);
        BoxData<double, DIM> x_i(b_i);
        BoxData<double, 1> J_i(b_i);
        map.apply(x_i, J_i, block);
        double J0 = J_i.absMax(); //J is a constant
        BoxData<double, C> phi_0 = forall_p<double, C>(f_phiM, block, x_i, k, offset);
        phi_i |= C2C(phi_0);

        BoxData<double, C> lphi_0 = forall_p<double, C>(f_LphiM, block, x_i, k, offset);
        auto& lphi_i = lphi[iter];
        lphi_i |= C2C(lphi_0);
        lphi_i *= J0;
    }
    auto& state = phi[*layout.begin()];
    auto& rhs = lphi[*layout.begin()];
    auto& X = map.map()[*layout.begin()];
    auto& J = map.jacobian()[*layout.begin()];
    auto dx = map.dx(0);
    
    std::cout << "X: " << X.box() << std::endl;
#if 1
    auto NT = Cofactor(X);
    std::cout << "NT: " << NT.box() << std::endl;

    auto M0 = Operator::cofactorMatrix(NT, 0);
    auto M1 = CofactorMatrix(NT, 0);

    M0.printData();
    M1.printData();
#else

    //------------------------------------------------------------------------------
    // Set Up Boundary Condition (Focusing on a single patch/block)
    Box B0 = layout[*layout.begin()];

    BCSet<double, C, HOST> bcs(B0);
    bcs.add(Face(0, Side::Lo), DIRICHLET);
    bcs.get(Face(0, Side::Lo)).data.setVal(0);
    bcs.add(Face(1, Side::Lo), DIRICHLET);
    bcs.get(Face(1, Side::Lo)).data.setVal(0);
    bcs.print();


    //------------------------------------------------------------------------------
    // Compute Gradient
    FluxBoxData<double, C, HOST> gradPhi(B0);
    Gradient(gradPhi, state, bcs, X, J, map.dx(0)); 
#endif
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
