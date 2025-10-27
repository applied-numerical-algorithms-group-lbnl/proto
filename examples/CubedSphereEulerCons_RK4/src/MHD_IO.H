#include "Proto.H"
#include "Inputs_Parsing.H"
#include "BoxOp_EulerCubedSphere.H"
#include "MHD_Constants.H"
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <cmath>
#include <limits>
#include <string>

time_t mktimeUTC(struct tm *time_tm)
{
  // *** enter in UTC mode
  char *oldTZ = getenv("TZ");
  putenv("TZ=UTC");
  tzset();
  // ***

  time_t ret = mktime(time_tm);

  // *** Restore previous TZ
  if (oldTZ == NULL)
  {
    putenv("TZ=");
  }
  else
  {
    char buff[255];
    sprintf(buff, "TZ=%s", oldTZ);
    putenv(buff);
  }
  tzset();
  // ***

  return ret;
}

time_t getPosixTime(double a_time)
{
  int year = (int)(floor(a_time));
  int daysYear = (year % 4 == 0 ? 366 : 365);

  tm time_tm;
  bzero(&time_tm, sizeof(tm));

  time_tm.tm_year = year - 1900;
  time_tm.tm_sec = (int)((a_time - year) * daysYear * 24.0 * 3600.0);
  time_tm.tm_mday = 1;
  time_tm.tm_isdst = 0;

  time_t posix_time = mktimeUTC(&time_tm);

  return posix_time;
}

double getPhysTime(double a_time)
{
  // double   timeRef      = eos_AU/m_lismV; // sec
  double BC_start_time = ParseInputs::get_BC_start_time();
  time_t startBC = getPosixTime(BC_start_time);
  double timeSec = (a_time); //*timeRef; // time in sec
  time_t timeSec_t = (time_t)(floor(timeSec));

  time_t curTime = startBC + timeSec_t;

  tm curTime_tm;
  // gmtime_r(&curTime,&curTime_tm);
  localtime_r(&curTime, &curTime_tm);

  int daysYear = (curTime_tm.tm_year % 4 == 0 ? 366 : 365);

  double PhysTime = curTime_tm.tm_year + 1900.0 +
                    (double)(curTime_tm.tm_yday) / (double)(daysYear) +
                    (double)(curTime_tm.tm_hour) / (double)(daysYear * 24.0) +
                    (double)(curTime_tm.tm_min) / (double)(daysYear * 24.0 * 60.0) +
                    (double)(curTime_tm.tm_sec + timeSec - timeSec_t) / (double)(daysYear * 24.0 * 3600.0);

  return PhysTime;
}

template<
  template<MemType> typename MAP,
  typename T, unsigned int C, MemType MEM, Centering CTR
>
void writeTecplotSlice_plane(
    const MBLevelMap<MAP, MEM>&            a_map,         // PR_NODE: x,y,z at nodes
    const MBLevelBoxData<T, C, MEM, CTR>&  a_data,        // cell-centered variables
    const std::vector<std::string>&        a_varNames,    // size C
    const std::string&                     a_filename,    // one Tecplot .dat
    double                                 a_physTime,
    double                                 a_time,
    int                                     a_iter, 
    int                                     a_axis,        // 0=x, 1=y, 2=z
    double                                  a_value = 0.0  // plane value (default 0)
)
{
    if (DIM != 3) {
        cout << "writePlaneSlice_Tecplot_SingleFile: DIM!=3 not supported\n";
        return;
    }
    if ((int)a_varNames.size() != (int)C) {
        cout << "writePlaneSlice_Tecplot_SingleFile: a_varNames.size()!=C\n";
        return;
    }
    if (a_axis < 0 || a_axis > 2) {
        cout << "writePlaneSlice_Tecplot_SingleFile: a_axis must be 0(x)/1(y)/2(z)\n";
        return;
    }

    const_cast<MBLevelMap<MAP, MEM>&>(a_map).initialize();

    std::ostringstream zbuf;
    zbuf.setf(std::ios::scientific);
    zbuf << std::setprecision(6);

    constexpr int Xc = 0, Yc = 1, Zc = 2;
    const int comp = (a_axis == 0 ? Xc : (a_axis == 1 ? Yc : Zc)); // component being sliced
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double tol = 1e-12;

    auto lerp = [](double a,double b,double t){ return a + t*(b-a); };

    auto titleAxis = [&](int ax)->const char*{
        return (ax==0 ? "x" : (ax==1 ? "y" : "z"));
    };

    int patchOrdinal = 0;
    for (auto it : a_data.layout()) {
        const auto& dataPatch = a_data[it];      // cell-centered
        const auto& mapPatch  = a_map.map()[it]; // node-centered

        // Cell interior box and sizes
        const Box  cellBox = dataPatch.box();
        const auto lo = cellBox.low();
        const auto hi = cellBox.high();
        const int ni = hi[0]-lo[0]+1;  // cells in i
        const int nj = hi[1]-lo[1]+1;  // cells in j
        const int nk = hi[2]-lo[2]+1;  // cells in k
        if (ni<=0 || nj<=0 || nk<=0) { ++patchOrdinal; continue; }

        // Skip patches with no plane crossing in their node envelope for chosen axis/value
        double vmin = +1e300, vmax = -1e300;
        {
            const Box mapBox = mapPatch.box();
            for (int k=lo[2]; k<=hi[2]+1; ++k)
            for (int j=lo[1]; j<=hi[1]+1; ++j)
            for (int i=lo[0]; i<=hi[0]+1; ++i) {
                Point pn(i,j,k);
                if (!cellBox.containsPoint(pn)) continue;
                const double v = mapPatch(pn, comp);
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
        }
        if (vmin > a_value + tol || vmax < a_value - tol) { ++patchOrdinal; continue; }

        const int blockId = a_data.layout().block(it);
        const bool isZface = (blockId == 0 || blockId == 1); // your mapping: 0/1 are ±Z

        // Choose marching direction and node grid based on your rules
        enum MarchDir { MARCH_I, MARCH_J, MARCH_K };
        MarchDir march;
        // Default grid sizes and lin indexers (set inside branches)
        int Nin=0, Njn=0, Nkn=0;

        // axis rules:
        // - Z plane: always march K (original behavior)
        // - Y plane: march J on non-±Z, march K on ±Z
        // - X plane: march I on non-±Z, march K on ±Z (rule-of-thumb)
        if (a_axis == Zc) {
            march = MARCH_K;
        } else if (a_axis == Yc) {
            march = isZface ? MARCH_K : MARCH_J;
        } else { // a_axis == Xc
            march = MARCH_J;
        }

        const Box  boxMap  = mapPatch.box();
        const auto mhi = boxMap.high();

        auto nodeClamped = [&](int ii, int jj, int kk) {
            int ci = std::min(ii, mhi[0]);
            int cj = std::min(jj, mhi[1]);
            int ck = std::min(kk, mhi[2]);
            return Point(ci, cj, ck);
        };

        // Common helpers for cell-center averages of map coords
        auto avg8_at = [&](int i,int j,int k,int component)->double{
            Point b00 = nodeClamped(lo[0]+i  , lo[1]+j  , lo[2]+k);
            Point b10 = nodeClamped(lo[0]+i+1, lo[1]+j  , lo[2]+k);
            Point b01 = nodeClamped(lo[0]+i  , lo[1]+j+1, lo[2]+k);
            Point b11 = nodeClamped(lo[0]+i+1, lo[1]+j+1, lo[2]+k);
            Point t00 = nodeClamped(lo[0]+i  , lo[1]+j  , lo[2]+k+1);
            Point t10 = nodeClamped(lo[0]+i+1, lo[1]+j  , lo[2]+k+1);
            Point t01 = nodeClamped(lo[0]+i  , lo[1]+j+1, lo[2]+k+1);
            Point t11 = nodeClamped(lo[0]+i+1, lo[1]+j+1, lo[2]+k+1);
            return 0.125*( mapPatch(b00,component)+mapPatch(b10,component)+mapPatch(b01,component)+mapPatch(b11,component)
                         + mapPatch(t00,component)+mapPatch(t10,component)+mapPatch(t01,component)+mapPatch(t11,component) );
        };

        // ============ Three branches: march K / march J / march I ============
        if (march == MARCH_K) {
            // Grid (i,j), cells ni*nj
            Nin = ni + 1; Njn = nj + 1;
            auto linN = [&](int i,int j){ return j*Nin + i; }; // nodes (i,j)
            auto linC = [&](int i,int j){ return j*ni  + i;  }; // cells (i,j)

            std::vector<double> X(Nin*Njn, NaN), Y(Nin*Njn, NaN), Z(Nin*Njn, NaN);
            // Intersect vertical node columns (vary k) with comp = a_value
            for (int jn=0; jn<Njn; ++jn) {
                for (int in=0; in<Nin; ++in) {
                    const int gi = lo[0] + in;
                    const int gj = lo[1] + jn;
                    const int mi = (gi > mhi[0]) ? mhi[0] : gi;
                    const int mj = (gj > mhi[1]) ? mhi[1] : gj;

                    bool hit=false;
                    for (int k=0; k<nk; ++k) {
                        const int kk  = lo[2] + k;
                        const int kk1 = lo[2] + std::min(k+1, nk);
                        Point p0 = nodeClamped(mi, mj, kk );
                        Point p1 = nodeClamped(mi, mj, kk1);
                        if (!boxMap.containsPoint(p0) || !boxMap.containsPoint(p1)) continue;

                        const double v0 = mapPatch(p0, comp);
                        const double v1 = mapPatch(p1, comp);
                        const bool cross = ( (v0 - a_value)<=+tol && (v1 - a_value)>=-tol ) ||
                                           ( (v0 - a_value)>=-tol && (v1 - a_value)<=+tol );
                        if (!cross) continue;
                        const double denom = (v1 - v0);
                        const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_value - v0)/denom);

                        const double x0 = mapPatch(p0, Xc), x1 = mapPatch(p1, Xc);
                        const double y0 = mapPatch(p0, Yc), y1 = mapPatch(p1, Yc);
                        const double z0 = mapPatch(p0, Zc), z1 = mapPatch(p1, Zc);

                        const int L = linN(in, jn);
                        X[L] = lerp(x0, x1, t);
                        Y[L] = lerp(y0, y1, t);
                        Z[L] = lerp(z0, z1, t);
                        // Pin sliced component numerically to a_value to avoid 1e-15 drift
                        if (comp==Xc) X[L]=a_value; else if (comp==Yc) Y[L]=a_value; else Z[L]=a_value;

                        hit = true; break;
                    }
                    (void)hit;
                }
            }

            // Cell-center interpolation along k
            std::vector<std::vector<double>> VarsCell(C, std::vector<double>(ni*nj, NaN));
            for (int j=0; j<nj; ++j) {
                for (int i=0; i<ni; ++i) {
                    bool have=false;
                    for (int k=0; k<nk-1; ++k) {
                        double c0 = avg8_at(i,j,k  , comp);
                        double c1 = avg8_at(i,j,k+1, comp);
                        const bool cross = ( (c0 - a_value)<=+tol && (c1 - a_value)>=-tol ) ||
                                           ( (c0 - a_value)>=-tol && (c1 - a_value)<=+tol );
                        if (!cross) continue;
                        const double denom = (c1 - c0);
                        const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_value - c0)/denom);

                        Point cp0(lo[0]+i, lo[1]+j, lo[2]+k  );
                        Point cp1(lo[0]+i, lo[1]+j, lo[2]+k+1);
                        for (unsigned c=0; c<C; ++c) {
                            const double v0 = static_cast<double>(dataPatch(cp0, c));
                            const double v1 = static_cast<double>(dataPatch(cp1, c));
                            VarsCell[c][linC(i,j)] = lerp(v0, v1, t);
                        }
                        have = true; break;
                    }
                    (void)have;
                }
            }

            // Emit zone (i,j)
            zbuf << "ZONE T=\"B" << a_data.layout().block(it) << "_P" << patchOrdinal
                 << "\" DATAPACKING=BLOCK, VARLOCATION=([4-" << (3+(int)C)
                 << "]=CELLCENTERED), I=" << (ni+1) << ", J=" << (nj+1) << "\n";

            auto dumpBlockN = [&](const std::vector<double>& A){
                for (int jj=0; jj<Njn; ++jj) {
                    for (int ii=0; ii<Nin; ++ii) zbuf << A[jj*Nin + ii] << " ";
                    zbuf << "\n";
                }
            };
            auto dumpBlockC = [&](const std::vector<double>& A){
                for (int jj=0; jj<nj; ++jj) {
                    for (int ii=0; ii<ni; ++ii) zbuf << A[jj*ni + ii] << " ";
                    zbuf << "\n";
                }
            };

            dumpBlockN(X); dumpBlockN(Y); dumpBlockN(Z);
            for (unsigned c=0; c<C; ++c) dumpBlockC(VarsCell[c]);
        }
        else if (march == MARCH_J) {
            // Grid (i,k), cells ni*nk
            Nin = ni + 1; Nkn = nk + 1;
            auto linN = [&](int i,int k){ return k*Nin + i; }; // nodes (i,k)
            auto linC = [&](int i,int k){ return k*ni  + i;  }; // cells (i,k)

            std::vector<double> X(Nin*Nkn, NaN), Y(Nin*Nkn, NaN), Z(Nin*Nkn, NaN);
            // Intersect along j at fixed (i,k): comp = a_value
            for (int kn=0; kn<Nkn; ++kn) {
                for (int in=0; in<Nin; ++in) {
                    const int gi = lo[0] + in;
                    const int gk = lo[2] + kn;
                    const int mi = (gi > mhi[0]) ? mhi[0] : gi;
                    const int mk = (gk > mhi[2]) ? mhi[2] : gk;

                    bool hit=false;
                    for (int j=0; j<nj; ++j) {
                        const int jj  = lo[1] + j;
                        const int jj1 = lo[1] + std::min(j+1, nj);
                        Point p0 = nodeClamped(mi, jj , mk);
                        Point p1 = nodeClamped(mi, jj1, mk);
                        if (!boxMap.containsPoint(p0) || !boxMap.containsPoint(p1)) continue;

                        const double v0 = mapPatch(p0, comp);
                        const double v1 = mapPatch(p1, comp);
                        const bool cross = ( (v0 - a_value)<=+tol && (v1 - a_value)>=-tol ) ||
                                           ( (v0 - a_value)>=-tol && (v1 - a_value)<=+tol );
                        if (!cross) continue;
                        const double denom = (v1 - v0);
                        const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_value - v0)/denom);

                        const double x0 = mapPatch(p0, Xc), x1 = mapPatch(p1, Xc);
                        const double y0 = mapPatch(p0, Yc), y1 = mapPatch(p1, Yc);
                        const double z0 = mapPatch(p0, Zc), z1 = mapPatch(p1, Zc);

                        const int L = linN(in, kn);
                        X[L] = lerp(x0, x1, t);
                        Y[L] = lerp(y0, y1, t);
                        Z[L] = lerp(z0, z1, t);
                        if (comp==Xc) X[L]=a_value; else if (comp==Yc) Y[L]=a_value; else Z[L]=a_value;

                        hit = true; break;
                    }
                    (void)hit;
                }
            }

            // Cell-center interpolation along j
            std::vector<std::vector<double>> VarsCell(C, std::vector<double>(ni*nk, NaN));
            for (int k=0; k<nk; ++k) {
                for (int i=0; i<ni; ++i) {
                    bool have=false;
                    for (int j=0; j<nj-1; ++j) {
                        double c0 = avg8_at(i,j  ,k, comp);
                        double c1 = avg8_at(i,j+1,k, comp);
                        const bool cross = ( (c0 - a_value)<=+tol && (c1 - a_value)>=-tol ) ||
                                           ( (c0 - a_value)>=-tol && (c1 - a_value)<=+tol );
                        if (!cross) continue;
                        const double denom = (c1 - c0);
                        const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_value - c0)/denom);

                        Point cp0(lo[0]+i, lo[1]+j  , lo[2]+k);
                        Point cp1(lo[0]+i, lo[1]+j+1, lo[2]+k);
                        for (unsigned c=0; c<C; ++c) {
                            const double v0 = static_cast<double>(dataPatch(cp0, c));
                            const double v1 = static_cast<double>(dataPatch(cp1, c));
                            VarsCell[c][linC(i,k)] = lerp(v0, v1, t);
                        }
                        have = true; break;
                    }
                    (void)have;
                }
            }

            // Emit zone (i,k)
            zbuf << "ZONE T=\"B" << a_data.layout().block(it) << "_P" << patchOrdinal
                 << "\" DATAPACKING=BLOCK, VARLOCATION=([4-" << (3+(int)C)
                 << "]=CELLCENTERED), I=" << (ni+1) << ", J=" << (nk+1) << "\n";

            auto dumpBlockN = [&](const std::vector<double>& A){
                for (int kk=0; kk<Nkn; ++kk) {
                    for (int ii=0; ii<Nin; ++ii) zbuf << A[kk*Nin + ii] << " ";
                    zbuf << "\n";
                }
            };
            auto dumpBlockC = [&](const std::vector<double>& A){
                for (int kk=0; kk<nk; ++kk) {
                    for (int ii=0; ii<ni; ++ii) zbuf << A[kk*ni + ii] << " ";
                    zbuf << "\n";
                }
            };

            dumpBlockN(X); dumpBlockN(Y); dumpBlockN(Z);
            for (unsigned c=0; c<C; ++c) dumpBlockC(VarsCell[c]);
        }
        else { // march == MARCH_I
            // Grid (j,k), cells nj*nk
            Njn = nj + 1; Nkn = nk + 1;
            auto linN = [&](int j,int k){ return k*Njn + j; }; // nodes (j,k)
            auto linC = [&](int j,int k){ return k*nj  + j;  }; // cells (j,k)

            std::vector<double> X(Njn*Nkn, NaN), Y(Njn*Nkn, NaN), Z(Njn*Nkn, NaN);
            // Intersect along i at fixed (j,k): comp = a_value
            for (int kn=0; kn<Nkn; ++kn) {
                for (int jn=0; jn<Njn; ++jn) {
                    const int gj = lo[1] + jn;
                    const int gk = lo[2] + kn;
                    const int mj = (gj > mhi[1]) ? mhi[1] : gj;
                    const int mk = (gk > mhi[2]) ? mhi[2] : gk;

                    bool hit=false;
                    for (int i=0; i<ni; ++i) {
                        const int ii  = lo[0] + i;
                        const int ii1 = lo[0] + std::min(i+1, ni);
                        Point p0 = nodeClamped(ii , mj, mk);
                        Point p1 = nodeClamped(ii1, mj, mk);
                        if (!boxMap.containsPoint(p0) || !boxMap.containsPoint(p1)) continue;

                        const double v0 = mapPatch(p0, comp);
                        const double v1 = mapPatch(p1, comp);
                        const bool cross = ( (v0 - a_value)<=+tol && (v1 - a_value)>=-tol ) ||
                                           ( (v0 - a_value)>=-tol && (v1 - a_value)<=+tol );
                        if (!cross) continue;
                        const double denom = (v1 - v0);
                        const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_value - v0)/denom);

                        const double x0 = mapPatch(p0, Xc), x1 = mapPatch(p1, Xc);
                        const double y0 = mapPatch(p0, Yc), y1 = mapPatch(p1, Yc);
                        const double z0 = mapPatch(p0, Zc), z1 = mapPatch(p1, Zc);

                        const int L = linN(jn, kn);
                        X[L] = lerp(x0, x1, t);
                        Y[L] = lerp(y0, y1, t);
                        Z[L] = lerp(z0, z1, t);
                        if (comp==Xc) X[L]=a_value; else if (comp==Yc) Y[L]=a_value; else Z[L]=a_value;

                        hit = true; break;
                    }
                    (void)hit;
                }
            }

            // Cell-center interpolation along i
            std::vector<std::vector<double>> VarsCell(C, std::vector<double>(nj*nk, NaN));
            for (int k=0; k<nk; ++k) {
                for (int j=0; j<nj; ++j) {
                    bool have=false;
                    for (int i=0; i<ni-1; ++i) {
                        double c0 = avg8_at(i  ,j,k, comp);
                        double c1 = avg8_at(i+1,j,k, comp);
                        const bool cross = ( (c0 - a_value)<=+tol && (c1 - a_value)>=-tol ) ||
                                           ( (c0 - a_value)>=-tol && (c1 - a_value)<=+tol );
                        if (!cross) continue;
                        const double denom = (c1 - c0);
                        const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_value - c0)/denom);

                        Point cp0(lo[0]+i  , lo[1]+j, lo[2]+k);
                        Point cp1(lo[0]+i+1, lo[1]+j, lo[2]+k);
                        for (unsigned c=0; c<C; ++c) {
                            const double v0 = static_cast<double>(dataPatch(cp0, c));
                            const double v1 = static_cast<double>(dataPatch(cp1, c));
                            VarsCell[c][linC(j,k)] = lerp(v0, v1, t);
                        }
                        have = true; break;
                    }
                    (void)have;
                }
            }

            // Emit zone (j,k)
            zbuf << "ZONE T=\"B" << a_data.layout().block(it) << "_P" << patchOrdinal
                 << "\" DATAPACKING=BLOCK, VARLOCATION=([4-" << (3+(int)C)
                 << "]=CELLCENTERED), I=" << (nj+1) << ", J=" << (nk+1) << "\n";

            auto dumpBlockN = [&](const std::vector<double>& A){
                for (int kk=0; kk<Nkn; ++kk) {
                    for (int jj=0; jj<Njn; ++jj) zbuf << A[kk*Njn + jj] << " ";
                    zbuf << "\n";
                }
            };
            auto dumpBlockC = [&](const std::vector<double>& A){
                for (int kk=0; kk<nk; ++kk) {
                    for (int jj=0; jj<nj; ++jj) zbuf << A[kk*nj + jj] << " ";
                    zbuf << "\n";
                }
            };

            dumpBlockN(X); dumpBlockN(Y); dumpBlockN(Z);
            for (unsigned c=0; c<C; ++c) dumpBlockC(VarsCell[c]);
        }

        ++patchOrdinal;
    }

    // ---- Gather and write single file on rank 0 ----
    std::string localStr = zbuf.str();
    std::vector<char> localBuf(localStr.begin(), localStr.end());
    int localN = (int)localBuf.size();

    int rank = 0, nprocs = 1;
#ifdef PR_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<int> counts, displs;
    if (rank == 0) counts.resize(nprocs);
    MPI_Gather(&localN, 1, MPI_INT, (rank==0?counts.data():nullptr), 1, MPI_INT, 0, comm);

    std::vector<char> big;
    if (rank == 0) {
        displs.resize(nprocs,0);
        int tot=0; for (int p=0;p<nprocs;++p){ displs[p]=tot; tot+=counts[p]; }
        big.resize(tot);
        MPI_Gatherv(localBuf.data(), localN, MPI_CHAR, big.data(), counts.data(), displs.data(), MPI_CHAR, 0, comm);

        std::ofstream ofs(a_filename, std::ios::out);
        ofs.setf(std::ios::scientific);
        ofs << std::setprecision(6);
        ofs << "TITLE = \"t=" << a_time << " step="<< a_iter<<"\"\n";
        ofs << "VARIABLES=\"X\" \"Y\" \"Z\"";
        for (const auto& nm : a_varNames) ofs << " \"" << nm << "\"";
        ofs << "\n";
        ofs << "DATASETAUXDATA PhysTime = \"" << std::setprecision(12) << a_physTime << "\"\n";
        ofs << std::setprecision(6);
        ofs.write(big.data(), (std::streamsize)big.size());
        ofs.close();
    } else {
        MPI_Gatherv(localBuf.data(), localN, MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR, 0, comm);
    }
#else
    std::ofstream ofs(a_filename, std::ios::out);
    ofs.setf(std::ios::scientific);
    ofs << std::setprecision(6);
    ofs << "TITLE = \"t=" << a_time << " step="<< a_iter<<"\"\n";
    ofs << "VARIABLES=\"X\" \"Y\" \"Z\"";
    for (const auto& nm : a_varNames) ofs << " \"" << nm << "\"";
    ofs << "\n";
    ofs << "DATASETAUXDATA PhysTime = \"" << std::setprecision(12) << a_physTime << "\"\n";
    ofs << std::setprecision(6);
    ofs << localStr;
    ofs.close();
#endif
}



// Extract a spherical slice r = a_radius into a single Tecplot file.
// Assumes "i" index marches monotonically outward in radius (increasing r).
// Robustness: we *search along i* for every (j,k) line, so even if r is not strictly
// monotone everywhere, any local crossing between i and i+1 will still be detected.
//
// - Geometry (X,Y,Z) is node-based on a (j,k) structured grid: I = nj+1, J = nk+1
// - Variables are interpolated at cell centers between (i,j,k) and (i+1,j,k): size nj*nk
//
// Usage:
// writeTecplotSlice_spherical(map, data, varNames, "r1.50.dat", time, /*radius=*/1.50);
template<
  template<MemType> typename MAP,
  typename T, unsigned int C, MemType MEM, Centering CTR
>
void writeTecplotSlice_spherical(
    const MBLevelMap<MAP, MEM>&            a_map,        // node-centered geometry: X,Y,Z
    const MBLevelBoxData<T, C, MEM, CTR>&  a_data,       // cell-centered variables
    const std::vector<std::string>&        a_varNames,   // size C
    const std::string&                     a_filename,   // output Tecplot .dat (single file)
    double                                  a_physTime,
    double                                 a_time,
    int                                     a_iter, 
    double                                  a_radius      // target radius r
)
{
    static_assert(DIM==3, "Tecplot slicer assumes DIM==3.");
    if ((int)a_varNames.size() != (int)C) {
        pr_out() << "writeSphericalSlice_Tecplot_SingleFile: a_varNames.size()!=C\n";
        return;
    }

    const_cast<MBLevelMap<MAP, MEM>&>(a_map).initialize();

    std::ostringstream zbuf;
    zbuf.setf(std::ios::scientific);
    zbuf << std::setprecision(6);

    constexpr int Xc = 0, Yc = 1, Zc = 2;
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double tol = 1e-12;

    auto lerp = [](double a,double b,double t){ return a + t*(b-a); };
    auto hypot3 = [](double x,double y,double z){ return std::sqrt(x*x + y*y + z*z); };

    int patchOrdinal = 0;
    for (auto it : a_data.layout()) {
        const auto& dataPatch = a_data[it];      // cell-centered
        const auto& mapPatch  = a_map.map()[it]; // node-centered

        // Cell interior box and sizes
        const Box  cellBox = dataPatch.box();
        const auto lo = cellBox.low();
        const auto hi = cellBox.high();
        const int ni = hi[0]-lo[0]+1;  // cells in i (radial march)
        const int nj = hi[1]-lo[1]+1;  // cells in j
        const int nk = hi[2]-lo[2]+1;  // cells in k
        if (ni<=0 || nj<=0 || nk<=0) { ++patchOrdinal; continue; }

        // Quick reject: check node-envelope radii against target
        double rmin = +1e300, rmax = -1e300;
        {
            const Box mapBox = mapPatch.box();
            for (int k=lo[2]; k<=hi[2]+1; ++k)
            for (int j=lo[1]; j<=hi[1]+1; ++j)
            for (int i=lo[0]; i<=hi[0]+1; ++i) {
                Point pn(i,j,k);
                if (!cellBox.containsPoint(pn)) continue; // consider overlap region
                const double x = mapPatch(pn, Xc);
                const double y = mapPatch(pn, Yc);
                const double z = mapPatch(pn, Zc);
                const double r = hypot3(x,y,z);
                if (r < rmin) rmin = r;
                if (r > rmax) rmax = r;
            }
        }
        if (rmin > a_radius + tol || rmax < a_radius - tol) { ++patchOrdinal; continue; }

        const int blockId = a_data.layout().block(it);
        (void)blockId; // (not needed here, but kept for symmetry with plane slicer)

        // Node grid for Tecplot geometry: use (j,k) since we march along i (radial)
        const int Njn = nj + 1, Nkn = nk + 1;
        auto linN = [&](int j,int k){ return k*Njn + j; }; // nodes (j,k)
        auto linC = [&](int j,int k){ return k*nj  + j;  }; // cells (j,k)

        std::vector<double> X(Njn*Nkn, NaN), Y(Njn*Nkn, NaN), Z(Njn*Nkn, NaN);

        const Box  boxMap  = mapPatch.box();
        const auto mhi = boxMap.high();

        auto nodeClamped = [&](int ii, int jj, int kk) {
            int ci = std::min(ii, mhi[0]);
            int cj = std::min(jj, mhi[1]);
            int ck = std::min(kk, mhi[2]);
            return Point(ci, cj, ck);
        };

        // ---------- NODE GEOMETRY: intersect along i for each (j,k) node-line ----------
        for (int kn=0; kn<Nkn; ++kn) {
            for (int jn=0; jn<Njn; ++jn) {
                const int gj = lo[1] + jn;                 // node j (may be hi[1]+1)
                const int gk = lo[2] + kn;                 // node k (may be hi[2]+1)
                const int mj = (gj > mhi[1]) ? mhi[1] : gj;
                const int mk = (gk > mhi[2]) ? mhi[2] : gk;

                bool hit=false;
                for (int i=0; i<ni; ++i) { // segment between i and i+1
                    const int ii  = lo[0] + i;
                    const int ii1 = lo[0] + std::min(i+1, ni);
                    Point p0 = nodeClamped(ii , mj, mk);
                    Point p1 = nodeClamped(ii1, mj, mk);
                    if (!boxMap.containsPoint(p0) || !boxMap.containsPoint(p1)) continue;

                    const double x0 = mapPatch(p0, Xc), x1 = mapPatch(p1, Xc);
                    const double y0 = mapPatch(p0, Yc), y1 = mapPatch(p1, Yc);
                    const double z0 = mapPatch(p0, Zc), z1 = mapPatch(p1, Zc);
                    const double r0 = hypot3(x0,y0,z0);
                    const double r1 = hypot3(x1,y1,z1);

                    const bool cross = ( (r0 - a_radius)<=+tol && (r1 - a_radius)>=-tol ) ||
                                       ( (r0 - a_radius)>=-tol && (r1 - a_radius)<=+tol );
                    if (!cross) continue;

                    const double denom = (r1 - r0);
                    const double t = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_radius - r0)/denom);

                    const int L = linN(jn, kn);
                    X[L] = lerp(x0, x1, t);
                    Y[L] = lerp(y0, y1, t);
                    Z[L] = lerp(z0, z1, t);
                    // Numerically nudge exactly onto the sphere (optional):
                    // const double s = a_radius / hypot3(X[L],Y[L],Z[L]); X[L]*=s; Y[L]*=s; Z[L]*=s;

                    hit = true; break;
                }
                (void)hit; // stays NaN if truly no intersection on this line
            }
        }

        // ---------- CELL-CENTER VARIABLES: interpolate along i at fixed (j,k) ----------
        // Average node map to cell-center (i,j,k) using the 8 corners
        auto avg8_at = [&](int i,int j,int k,int comp)->double{
            Point b00 = nodeClamped(lo[0]+i  , lo[1]+j  , lo[2]+k);
            Point b10 = nodeClamped(lo[0]+i+1, lo[1]+j  , lo[2]+k);
            Point b01 = nodeClamped(lo[0]+i  , lo[1]+j+1, lo[2]+k);
            Point b11 = nodeClamped(lo[0]+i+1, lo[1]+j+1, lo[2]+k);
            Point t00 = nodeClamped(lo[0]+i  , lo[1]+j  , lo[2]+k+1);
            Point t10 = nodeClamped(lo[0]+i+1, lo[1]+j  , lo[2]+k+1);
            Point t01 = nodeClamped(lo[0]+i  , lo[1]+j+1, lo[2]+k+1);
            Point t11 = nodeClamped(lo[0]+i+1, lo[1]+j+1, lo[2]+k+1);
            return 0.125*( mapPatch(b00,comp)+mapPatch(b10,comp)+mapPatch(b01,comp)+mapPatch(b11,comp)
                         + mapPatch(t00,comp)+mapPatch(t10,comp)+mapPatch(t01,comp)+mapPatch(t11,comp) );
        };

        // Helper to get center radius at (i,j,k)
        auto center_r = [&](int i,int j,int k)->double{
            const double xc = avg8_at(i,j,k, Xc);
            const double yc = avg8_at(i,j,k, Yc);
            const double zc = avg8_at(i,j,k, Zc);
            return hypot3(xc,yc,zc);
        };

        std::vector<std::vector<double>> VarsCell(C, std::vector<double>(nj*nk, NaN));

        for (int k=0; k<nk; ++k) {
            for (int j=0; j<nj; ++j) {
                bool have = false;
                for (int i=0; i<ni-1; ++i) {
                    const double r0 = center_r(i  , j, k);
                    const double r1 = center_r(i+1, j, k);

                    const bool cross = ( (r0 - a_radius)<=+tol && (r1 - a_radius)>=-tol ) ||
                                       ( (r0 - a_radius)>=-tol && (r1 - a_radius)<=+tol );
                    if (!cross) continue;

                    const double denom = (r1 - r0);
                    const double tcc   = (std::fabs(denom) < 1e-30) ? 0.0 : ((a_radius - r0)/denom);

                    Point cp0(lo[0]+i  , lo[1]+j, lo[2]+k);
                    Point cp1(lo[0]+i+1, lo[1]+j, lo[2]+k);
                    for (unsigned c=0; c<C; ++c) {
                        const double v0 = static_cast<double>(dataPatch(cp0, c));
                        const double v1 = static_cast<double>(dataPatch(cp1, c));
                        VarsCell[c][linC(j,k)] = lerp(v0, v1, tcc);
                    }
                    have = true;
                    break;
                }
                (void)have; // remains NaN where sphere truly does not cut the column
            }
        }

        // ---- Emit ZONE for this patch (structured (j,k)) ----
        zbuf << "ZONE T=\"B" << a_data.layout().block(it) << "_P" << patchOrdinal
             << "\" DATAPACKING=BLOCK, VARLOCATION=([4-" << (3+(int)C)
             << "]=CELLCENTERED), I=" << (nj+1) << ", J=" << (nk+1) << "\n";

        // Dump helpers
        auto dumpBlockN = [&](const std::vector<double>& A){
            for (int kk=0; kk<Nkn; ++kk) {
                for (int jj=0; jj<Njn; ++jj) zbuf << A[linN(jj,kk)] << " ";
                zbuf << "\n";
            }
        };
        auto dumpBlockC = [&](const std::vector<double>& A){
            for (int kk=0; kk<nk; ++kk) {
                for (int jj=0; jj<nj; ++jj) zbuf << A[linC(jj,kk)] << " ";
                zbuf << "\n";
            }
        };

        dumpBlockN(X);
        dumpBlockN(Y);
        dumpBlockN(Z);
        for (unsigned c=0; c<C; ++c) dumpBlockC(VarsCell[c]);

        ++patchOrdinal;
    }

    // ---- Gather and write single file on rank 0 ----
    std::string localStr = zbuf.str();
    std::vector<char> localBuf(localStr.begin(), localStr.end());
    int localN = (int)localBuf.size();

    int rank = 0, nprocs = 1;
#ifdef PR_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<int> counts, displs;
    if (rank == 0) counts.resize(nprocs);
    MPI_Gather(&localN, 1, MPI_INT, (rank==0?counts.data():nullptr), 1, MPI_INT, 0, comm);

    std::vector<char> big;
    if (rank == 0) {
        displs.resize(nprocs,0);
        int tot=0; for (int p=0;p<nprocs;++p){ displs[p]=tot; tot+=counts[p]; }
        big.resize(tot);
        MPI_Gatherv(localBuf.data(), localN, MPI_CHAR, big.data(), counts.data(), displs.data(), MPI_CHAR, 0, comm);

        std::ofstream ofs(a_filename, std::ios::out);
        ofs.setf(std::ios::scientific);
        ofs << std::setprecision(6);
        ofs << "TITLE = \"t=" << a_time << " step="<< a_iter<<"\"\n";
        ofs << "VARIABLES=\"X\" \"Y\" \"Z\"";
        for (const auto& nm : a_varNames) ofs << " \"" << nm << "\"";
        ofs << "\n";
        ofs << "DATASETAUXDATA PhysTime = \"" << std::setprecision(12) << a_physTime << "\"\n";
        ofs << std::setprecision(6);
        ofs.write(big.data(), (std::streamsize)big.size());
        ofs.close();
    } else {
        MPI_Gatherv(localBuf.data(), localN, MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR, 0, comm);
    }
#else
    std::ofstream ofs(a_filename, std::ios::out);
    ofs.setf(std::ios::scientific);
    ofs << std::setprecision(6);
    ofs << "TITLE = \"t=" << a_time << " step="<< a_iter<<"\"\n";
    ofs << "VARIABLES=\"X\" \"Y\" \"Z\"";
    for (const auto& nm : a_varNames) ofs << " \"" << nm << "\"";
    ofs << "\n";
    ofs << "DATASETAUXDATA PhysTime = \"" << std::setprecision(12) << a_physTime << "\"\n";
    ofs << std::setprecision(6);
    ofs << localStr;
    ofs.close();
#endif
}

template<typename BOX_OP>
void Write_W(MBLevelBoxData<double, NUMCOMPS, HOST> &a_JU,
             BOX_OP &eulerOp,
             MBInterpOp &a_iop,
             int iter,
             double time = 0.0,
             double dt = 0.0,
             bool write_data = false,
             bool write_slice = false)
{
  double gamma = ParseInputs::get_gamma();
  HDF5Handler h5;
  a_JU.exchange();
  h5.setTime(time);
	h5.setTimestep(dt);
  h5.setIter(iter);
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  auto layout = a_JU.layout();
  auto map = CubedSphereShell::Map(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS+2, HOST> Wout(layout, Point::Zeros());
  MBLevelBoxData<double, NUMCOMPS, HOST> Uout(layout, Point::Zeros());
  MBLevelBoxData<double, NUMCOMPS, HOST> W_cart_out(layout, Point::Zeros());
  MBLevelBoxData<double, NUMCOMPS, HOST> W_sph_out(layout, Point::Zeros());
  MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
  MBLevelBoxData<double, 1, HOST> radiusLev(layout, OP::ghost() + 2);
  for (auto dit : layout)
  {
    BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
    BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
    eulerOp[dit].radialMetrics(radiusLev[dit], Dr, adjDr, dVolrLev[dit], Dr.box());
  }
  Point size = layout.domain().getBlock(0).sizes();
  double dxiPerp = 1.0/size[1];
  
  bool infOrNan = false;
  int scale_data_with_r = ParseInputs::get_scale_data_with_r(); 
  Array<std::string,NUMCOMPS+2> names = {"density", "VxPerp", "VyPerp", "VzPerp", "P", "BxPerp", "ByPerp", "BzPerp","Vr","Br*r^2"};//,"e^entropy"};
  if (scale_data_with_r == 1){
    names[0] = "density_scaled";
    names[4] = "P_scaled";
    names[5] = "BxPerp_scaled";
    names[6] = "ByPerp_scaled";
    names[7] = "BzPerp_scaled";
    names[9] = "Br_scaled"; 
  }
  Array<double, DIM> dx;
  for (auto dit : layout)
  {
    auto& JU_i = a_JU[dit];
    BoxData<double,1,HOST> deta = CubedSphereShell::detA<double,HOST>(JU_i.box(),dxiPerp,0); 
    BoxData<double,NUMCOMPS+2,HOST> USemiSph(layout[dit].grow(Point::Basis(0,2)));
    auto& Wout_i = Wout[dit];
    USemiSph.setToZero();
    unsigned int block = layout.block(dit);
    Box blockBox = layout.getBlock(block).domain().box();
    auto& dVolr_i = dVolrLev[dit];
    Box bx_i = layout[dit];
    CubedSphereShell::
      consToSemiSphNGEuler(USemiSph,JU_i,dVolr_i,bx_i,
                           blockBox,dxiPerp,block,4);
    dx = eulerOp[dit].dx();
    double half = 0.5;
    BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
            (f_cubedSphereMap3,radiusLev[dit].box(),radiusLev[dit],dx,half,half,block);  
    forallInPlace([ ] PROTO_LAMBDA
    (Var<double,NUMCOMPS+2,HOST>& a_Wout,
     Var<double,NUMCOMPS,HOST>& a_W_sph_out,
     const Var<double,NUMCOMPS+2,HOST>& a_USemiSph,
     double a_gamma,
     const Var<double,1,HOST>& a_radius,
     const Var<double,DIM,HOST>& a_Xcart,
     int a_scale_data_with_r)
    {
      // Primitive variables with vec_non_rad as vec_x_non, vec_y_non, vec_z_non
      a_Wout(0) = a_USemiSph(0);
      a_Wout(1) = a_USemiSph(1) / a_Wout(0);
      a_Wout(2) = a_USemiSph(2) / a_Wout(0);
      a_Wout(3) = a_USemiSph(3) / a_Wout(0);
      a_Wout(5) = a_USemiSph(5);
      a_Wout(6) = a_USemiSph(6);
      a_Wout(7) = a_USemiSph(7);
      a_Wout(8) = a_USemiSph(8) / a_Wout(0);
      a_Wout(9) = a_USemiSph(9);
      double ke = a_Wout(0)*(a_Wout(1)*a_Wout(1)
                         + a_Wout(2)*a_Wout(2)
                         + a_Wout(3)*a_Wout(3) + a_Wout(8)*a_Wout(8))*.5;
      double me = (a_Wout(5)*a_Wout(5)
              + a_Wout(6)*a_Wout(6)
              + a_Wout(7)*a_Wout(7) + a_Wout(9)*a_Wout(9))/(8.0*M_PI);
      a_Wout(4) = (a_USemiSph(4) - ke - me)*(a_gamma - 1.0);


      a_W_sph_out(iRHO) = a_Wout(0);
      a_W_sph_out(iP) = a_Wout(4);
      // Convert non rad cartesian to spherical velocities
      double x = a_Xcart(0);
      double y = a_Xcart(1);
      double z = a_Xcart(2);
      double r = sqrt(x*x + y*y + z*z);
      double theta = acos(z/r);
      double phi = atan2(y,x);
      // Convert Cartesian to spherical velocities
      a_W_sph_out(iVX) = a_Wout(8);
      a_W_sph_out(iVY) = cos(theta)*cos(phi)*a_Wout(1) + cos(theta)*sin(phi)*a_Wout(2) - sin(theta)*a_Wout(3);
      a_W_sph_out(iVZ) = -sin(phi)*a_Wout(1) + cos(phi)*a_Wout(2);
      // Convert Cartesian to spherical magnetic fields
      a_W_sph_out(iBX) = a_Wout(9);
      a_W_sph_out(iBY) = cos(theta)*cos(phi)*a_Wout(5) + cos(theta)*sin(phi)*a_Wout(6) - sin(theta)*a_Wout(7);
      a_W_sph_out(iBZ) = -sin(phi)*a_Wout(5) + cos(phi)*a_Wout(6);

      if (a_scale_data_with_r == 1){
        a_W_sph_out(iRHO) = a_W_sph_out(iRHO)*a_radius(0)*a_radius(0)/c_AU/c_AU/c_MP; // /cm^3 , scaled by r^2;
        a_W_sph_out(iVX) = a_W_sph_out(iVX)/1e5; // km/s;
        a_W_sph_out(iVY) = a_W_sph_out(iVY)/1e5; // km/s;
        a_W_sph_out(iVZ) = a_W_sph_out(iVZ)/1e5; // km/s;
        a_W_sph_out(iP) = a_W_sph_out(iP)*pow(a_radius(0)/c_AU,a_gamma*2)*1e12;    
        // scaled by r*(2*gamma)
        if (a_Wout(iP) < 0.) a_Wout(iP) = 0.;
        a_W_sph_out(iBX) = a_W_sph_out(iBX)*(a_radius(0)*a_radius(0)/c_AU/c_AU)*1e6; // microGauss
        a_W_sph_out(iBY) = a_W_sph_out(iBY)*(a_radius(0)/c_AU)*1e6; // microGauss
        a_W_sph_out(iBZ) = a_W_sph_out(iBZ)*(a_radius(0)/c_AU)*1e6; // microGauss
      }
      
    }, Wout[dit], W_sph_out[dit], USemiSph, gamma, radiusLev[dit], XCart, scale_data_with_r);

#if 0    
    forallInPlace([ ] PROTO_LAMBDA
      (Var<double,NUMCOMPS,HOST>& l_U,
        const Var<double,NUMCOMPS,HOST>& l_JU,
        const Var<double,1,HOST>& l_deta,
        const Var<double,1,HOST>& l_dVolr)
      {                    
        for (int comp = 0; comp < NUMCOMPS; comp++)
          {
            l_U(comp) = l_JU(comp)/(l_deta(0)*l_dVolr(0));
          } 
      },Uout[dit],JU_i,deta,dVolr_i);
#endif

#ifndef PR_DISABLE_INIT_CHECKS
    // bool infOrNanBlock = Wout[dit].containsInfOrNAN();
    // if (infOrNanBlock)
    //   {
    //     h5.writePatch(1.0,Wout[dit],"MHD_IO_InfOrNAN"+to_string(block));
    //     infOrNan = infOrNan || infOrNanBlock;
    //     cout << "Inf or Nan at Point " << Wout[dit].findInfOrNAN()
    //          << ", block = " << block << endl;
    //   }
#endif
#if 0
    if (procID() == 0)
    {
      Box bx = layout[dit];
      std::string output1D_file = ParseInputs::get_data_file_prefix()+"index" + "1D_" + std::to_string(iter)+".curve";
      std::ofstream ofs (output1D_file, std::ofstream::out);
      Point ptcen = bx.high()/2;
      int highxi0 = bx.high()[0];
      ptcen[0] = 0;
      
      for (int comp = 0; comp < NUMCOMPS+3; comp++)
        {
          ofs << "#" + names[comp] << endl;
          for (int ind = 0; ind < highxi0; ind++)
            {
              Point ptindex = ptcen;
              ptindex[0] = ind;
              Point ptindex_p = ptcen;
              ptindex_p[0] = ind+1;
              double r = (radiusLev[dit](ptindex,0) + radiusLev[dit](ptindex_p,0))*.5;
              if (comp < NUMCOMPS+2)
                {
                  double val = Wout[dit](ptindex,comp);
                  ofs << r << " " << val << endl;
                }
              else if (comp == NUMCOMPS+2)
                {
                  double p = Wout[dit](ptindex,4);
                  double rho = Wout[dit](ptindex,0);
                  double val = p/pow(rho,gamma);
                  ofs << r << " " << val << endl;
                }
            }
          ofs << endl;
        }
    }
#endif
  }
  
  if (write_data){
    std::string output_file = ParseInputs::get_data_file_prefix() + "_" + std::to_string(iter); 
    h5.writeMBLevel({"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"}, map, W_sph_out, output_file);
    if (procID() == 0) cout << "Data file written at iteration " << iter << endl;
    // Print this for more detailed output
    // h5.writeMBLevel({"density", "VxPerp", "VyPerp", "VzPerp", "P", "BxPerp", "ByPerp", "BzPerp","Vr","Br"}, map, Wout, "W_" + to_string(iter));
  }

  if (write_slice){
    double physical_time = getPhysTime(time);
    std::string tecplot_slice_file_x = ParseInputs::get_data_file_prefix() + ".x." + std::to_string(iter) + ".dat";
    std::string tecplot_slice_file_y = ParseInputs::get_data_file_prefix() + ".y." + std::to_string(iter) + ".dat";
    std::string tecplot_slice_file_z = ParseInputs::get_data_file_prefix() + ".z." + std::to_string(iter) + ".dat";
    std::string tecplot_spherical_slice_file = ParseInputs::get_data_file_prefix() + ".r." + std::to_string(iter) + ".dat";
    writeTecplotSlice_plane(map, W_sph_out, {"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"}, tecplot_slice_file_x, physical_time,time,iter, 0, 0.0);
    writeTecplotSlice_plane(map, W_sph_out, {"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"}, tecplot_slice_file_y, physical_time,time,iter, 1, 0.0);
    writeTecplotSlice_plane(map, W_sph_out, {"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"}, tecplot_slice_file_z, physical_time,time,iter, 2, 0.0);
    writeTecplotSlice_spherical(map, W_sph_out, {"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"}, tecplot_spherical_slice_file, physical_time, time, iter, /*radius=*/1.5e13); 
    if (procID() == 0) cout << "Tecplot slice files written at iteration " << iter << endl;                                     
  }
  
}

  void Write_Checkpoint(MBLevelBoxData<double, NUMCOMPS, HOST> &a_JU,
             int iter,
             int restart_step,
             double time = 0.0,
             double dt = 0.0)
{
  HDF5Handler h5;
  std::string check_file_name = ParseInputs::get_checkpoint_file_prefix() + "_" + std::to_string(iter);
  h5.setTime(time);
  h5.setTimestep(dt);
  h5.setIter(iter);
  h5.writeMBLevel(a_JU, check_file_name);
  if (procID() == 0) cout << "Checkpoint file written at iteration " << iter << endl;
  int checkpoint_cadence = ParseInputs::get_checkpoint_cadence();
  int iter_to_delete = iter - (checkpoint_cadence*ParseInputs::get_max_checkpoint_files());
  if (iter_to_delete > restart_step)
  {
    std::string filename_to_delete=ParseInputs::get_checkpoint_file_prefix()+"_"+std::to_string(iter_to_delete)+".hdf5";
    const char* str = filename_to_delete.c_str();
    if (procID() == 0) std::remove(str);
  }
}

// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
double interpolate(vector<double> &xData, vector<double> &yData, double x, bool extrapolate)
{
  int size = xData.size();

  int i = 0;                // find left end of interval for interpolation
  if (x >= xData[size - 2]) // special case: beyond right end
  {
    i = size - 2;
  }
  else
  {
    while (x > xData[i + 1])
      i++;
  }
  double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1]; // points on either side (unless beyond ends)
  if (!extrapolate)                                                          // if beyond ends of array and not extrapolating
  {
    if (x < xL)
      yR = yL;
    if (x > xR)
      yL = yR;
  }

  double dydx = (yR - yL) / (xR - xL); // gradient

  return yL + dydx * (x - xL); // linear interpolation
}



template<typename BOX_OP>
void Probe(const MBLevelBoxData<double, NUMCOMPS, HOST> &a_JU,
          const MBLevelMap<MBMap_CubedSphereShell, HOST>& a_map,
          BOX_OP &a_eulerOp,
          MBInterpOp &a_iop,
          int a_iter,
          double a_time,
          Array<double, DIM> a_dx,
          bool give_space = false)

{

  int pid = procID();
  ofstream outputFile;
  string Probe_data_file = ParseInputs::get_Probe_data_file();
  outputFile.open(Probe_data_file, std::ios::app);
  if (pid == 0 && give_space)
    outputFile << endl;
  double gamma = ParseInputs::get_gamma();
  std::string probe_file = ParseInputs::get_Probe_trajectory_file();

  int number_of_lines = 0;
  int totvars = 5; // Year DOY Radius Latitude Longitude
  std::string line;
  std::ifstream myfile(probe_file);

  PROTO_ASSERT(myfile.is_open(), "Error: Unable to open trajectory file");

  while (std::getline(myfile, line))
    ++number_of_lines;

  int row, col;
  double *my_array;
  my_array = new double[number_of_lines * totvars];

  ifstream pFile(probe_file);
  row = 0;
  while (!pFile.eof())
  {
    getline(pFile, line);
    stringstream ss(line);
    col = 0;

    while (ss >> my_array[row * totvars + col])
    {

      col++;
    }
    row++;
  }
  pFile.close();

  vector<double> time_traj_file, r_traj_file, lat_traj_file, lon_traj_file;
  int year, day;
  double rad_AU, latitude, longitude, relative_time;
  bool leap_year;
  for (int i = 0; i < number_of_lines; i++)
  {
    year = my_array[i * totvars + 0];
    day = my_array[i * totvars + 1];
    rad_AU = my_array[i * totvars + 2];
    latitude = my_array[i * totvars + 3];
    longitude = my_array[i * totvars + 4];
    leap_year = false;
    if ((year % 4) == 0)
    {
      leap_year = true;
      if ((year % 100) == 0 && (year % 400) != 0)
        leap_year = false;
    }
    if (leap_year)
    {
      // need to check with Tae if the values are at start of day or end of day. Will use (day-1) then.
      relative_time = year + (day - 1) / 366.0;
    }
    else
    {
      relative_time = year + (day - 1) / 365.0;
    }

    time_traj_file.push_back(relative_time);
    r_traj_file.push_back(rad_AU);
    lat_traj_file.push_back(latitude);
    lon_traj_file.push_back(longitude);
  }
  delete[] my_array;

  double physical_time = getPhysTime(a_time);
  double r_now = interpolate(time_traj_file, r_traj_file, physical_time, false);
  double latitude_now = interpolate(time_traj_file, lat_traj_file, physical_time, false);
  double longitude_now = interpolate(time_traj_file, lon_traj_file, physical_time, false);
  r_now *= c_AU;
  latitude_now = (90.0 - latitude_now);   // To co-latitude
  longitude_now = longitude_now + 181.02; // HGI to Carrington
  if (longitude_now >= 360.0)
    longitude_now -= 360;
  // longitude_now = longitude_now - 1.792951949248379e+02; // HGI to Carrington
  if (longitude_now < 0.0)
    longitude_now += 360;

  double x_n = sin(latitude_now * M_PI / 180) * cos(longitude_now * M_PI / 180);
  double y_n = sin(latitude_now * M_PI / 180) * sin(longitude_now * M_PI / 180);
  double z_n = cos(latitude_now * M_PI / 180);

  double x = r_now * x_n;
  double y = r_now * y_n;
  double z = r_now * z_n;

  double H3_to_RTN_00 = x_n;
  double H3_to_RTN_01 = y_n;
  double H3_to_RTN_02 = z_n;
  double H3_to_RTN_10 = -y_n;
  double H3_to_RTN_11 = x_n;
  double H3_to_RTN_12 = 0.0;
  double H3_to_RTN_20 = -x_n * z_n;
  double H3_to_RTN_21 = -y_n * z_n;
  double H3_to_RTN_22 = x_n * x_n + y_n * y_n;

  const Array<double, DIM> Xcart = {x, y, z};
  Array<double, DIM> Xi;
  // Array<double, DIM> Xi = a_map.inverse(Xcart); This inverse will work in future

  const auto& layout = a_map.layout();
  for (auto iter : layout)
  {
    BlockIndex block = layout.block(iter);
    auto& mapOp = a_map.op(block);
    Xi = mapOp.inverse(Xcart);
  }

  double eta0 = Xi[0];
  double eta1 = Xi[1];
  double eta2 = Xi[2];

  double pt0 = eta0 / a_dx[0];
  double pt1 = eta1 / a_dx[1];
  double pt2 = eta2 / a_dx[2];

  int pt0_nearest, pt1_nearest, pt2_nearest;
  int pt0_neigbor, pt1_neigbor, pt2_neigbor;

  pt0_nearest = floor(pt0);
  pt1_nearest = floor(pt1);
  pt2_nearest = floor(pt2);

  
  if (round(pt0) == floor(pt0))
  {
    pt0_neigbor = pt0_nearest - 1;
  }
  else
  {
    pt0_neigbor = pt0_nearest + 1;
  }

  if (round(pt1) == floor(pt1))
  {
    pt1_neigbor = pt1_nearest - 1;
  }
  else
  {
    pt1_neigbor = pt1_nearest + 1;
  }

  if (round(pt2) == floor(pt2))
  {
    pt2_neigbor = pt2_nearest - 1;
  }
  else
  {
    pt2_neigbor = pt2_nearest + 1;
  }

  double probed_values[NUMCOMPS];
  double probed_values_primitive[NUMCOMPS];
  Point index_cc(pt0_nearest, pt1_nearest, pt2_nearest);
  Point index_n0(pt0_neigbor, pt1_nearest, pt2_nearest);
  Point index_n1(pt0_nearest, pt1_neigbor, pt2_nearest);
  Point index_n2(pt0_nearest, pt1_nearest, pt2_neigbor);


  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;

  MBLevelBoxData<double, NUMCOMPS, HOST> JUTemp(layout, OP::ghost());
  MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
  MBLevelBoxData<double, 1, HOST> radiusLev(layout, OP::ghost() + 2);
  for (auto dit : a_JU.layout())
  {
    BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
    BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
    a_eulerOp[dit].radialMetrics(radiusLev[dit], Dr, adjDr, dVolrLev[dit], Dr.box());
    a_JU[dit].copyTo(JUTemp[dit]);
  }

  CubedSphereShell::consToCartInterpEuler(JUTemp, a_iop, dVolrLev, 4);


  
  for (auto dit : layout)
  {
    auto block = layout.block(dit);
    Box dbx0 = a_JU[dit].box();
    Box dbx1 = dbx0.grow(-NGHOST);
    BoxData<double, NUMCOMPS, HOST> U_dim(dbx0);
    JUTemp[dit].copyTo(U_dim);

    if (dbx1.containsPoint(index_cc))
    {
      for (int i = 0; i < NUMCOMPS; i++)
      {
        probed_values[i] = U_dim(index_cc, i);
        if (dbx0.containsPoint(index_n0))
          probed_values[i] += (U_dim(index_n0, i) - U_dim(index_cc, i)) * (pt0 - (pt0_nearest + 0.5)) / (pt0_neigbor - (pt0_nearest + 0.5));
        if (dbx0.containsPoint(index_n1))
          probed_values[i] += (U_dim(index_n1, i) - U_dim(index_cc, i)) * (pt1 - (pt1_nearest + 0.5)) / (pt1_neigbor - (pt1_nearest + 0.5));
        if (dbx0.containsPoint(index_n2))
          probed_values[i] += (U_dim(index_n2, i) - U_dim(index_cc, i)) * (pt2 - (pt2_nearest + 0.5)) / (pt2_neigbor - (pt2_nearest + 0.5));
      }

      double v2 = 0.0;
      double B2 = 0.0;
      probed_values_primitive[0] = probed_values[0];

      for (int i = 1; i <= DIM; i++)
      {
        double v, B;
        v = probed_values[i] / probed_values[0];
        B = probed_values[DIM + 1 + i];
        probed_values_primitive[i] = v;
        probed_values_primitive[DIM + 1 + i] = B;
        v2 += v * v;
        B2 += B * B;
      }

      probed_values_primitive[DIM + 1] = (probed_values[DIM + 1] - .5 * probed_values[0] * v2 - B2 / 8.0 / M_PI) * (gamma - 1.0);

      double rho = probed_values_primitive[0] / c_MP; // /cm^3;
      double Vx = probed_values_primitive[1] / 1e5;   // km/s;
      double Vy = probed_values_primitive[2] / 1e5;   // km/s;
      double Vz = probed_values_primitive[3] / 1e5;   // km/s;
      double p = probed_values_primitive[4] * 1e12;   // picodyne/cm*2
      double Bx = probed_values_primitive[5] * 1e6;   // microGauss
      double By = probed_values_primitive[6] * 1e6;   // microGauss
      double Bz = probed_values_primitive[7] * 1e6;   // microGauss

      double VR, VT, VN, BR, BT, BN;

      VR = H3_to_RTN_00 * Vx + H3_to_RTN_01 * Vy + H3_to_RTN_02 * Vz;
      VT = H3_to_RTN_10 * Vx + H3_to_RTN_11 * Vy + H3_to_RTN_12 * Vz;
      VN = H3_to_RTN_20 * Vx + H3_to_RTN_21 * Vy + H3_to_RTN_22 * Vz;

      BR = H3_to_RTN_00 * Bx + H3_to_RTN_01 * By + H3_to_RTN_02 * Bz;
      BT = H3_to_RTN_10 * Bx + H3_to_RTN_11 * By + H3_to_RTN_12 * Bz;
      BN = H3_to_RTN_20 * Bx + H3_to_RTN_21 * By + H3_to_RTN_22 * Bz;

      // if (probed_values[0] != 0) {
      outputFile << setw(16) << setprecision(12) << physical_time
                 << setw(11) << setprecision(4) << r_now / c_AU
                 << setw(11) << setprecision(4) << 90 - latitude_now
                 << setw(11) << setprecision(4) << longitude_now
                 << setw(11) << setprecision(4) << rho
                 << setw(11) << setprecision(4) << VR
                 << setw(11) << setprecision(4) << VT
                 << setw(11) << setprecision(4) << VN
                 << setw(11) << setprecision(4) << p
                 << setw(11) << setprecision(4) << BR
                 << setw(11) << setprecision(4) << BT
                 << setw(11) << setprecision(4) << BN
                 << endl;
      // }
    }
  }

  if (procID() == 0) cout << "Probed at iteration " << a_iter << endl;

  outputFile.close();
}
