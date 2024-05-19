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

template<typename BOX_OP>
void Write_W(MBLevelBoxData<double, NUMCOMPS, HOST> &a_JU,
             BOX_OP &eulerOp,
             MBInterpOp &a_iop,
             int iter,
             double time = 0.0,
             double dt = 0.0)
{
  double gamma = ParseInputs::get_gamma();
  HDF5Handler h5;
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  auto layout = a_JU.layout();
  auto map = CubedSphereShell::Map(layout, OP::ghost());

  MBLevelBoxData<double, NUMCOMPS, HOST> USph(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> JUTemp(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> Wout(layout, OP::ghost());
  MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
  MBLevelBoxData<double, 1, HOST> radiusLev(layout, OP::ghost() + 2);
  for (auto dit : a_JU.layout())
  {
    BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
    BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
    eulerOp[dit].radialMetrics(radiusLev[dit], Dr, adjDr, dVolrLev[dit], Dr.box());
    a_JU[dit].copyTo(JUTemp[dit]);
  }

  CubedSphereShell::consToSphInterpEuler(JUTemp, a_iop, dVolrLev, 4);

  for (auto dit : a_JU.layout())
  {
    State W;
    // State WBar;
    // eulerOp[dit].sphToPrim(W, WBar, JUTemp[dit], layout.block(dit));
    // Second order for output purposes. This then doesn't need BCs.
    eulerOp[dit].sphToPrim_2O(W, JUTemp[dit], layout.block(dit));
    forallInPlace([ ] PROTO_LAMBDA
    (Var<double,NUMCOMPS,HOST>& a_Wout,
    const Var<double,NUMCOMPS,HOST>& a_W)
    {
      a_Wout(0) = a_W(0)/c_MP; // /cm^3;
      a_Wout(1) = a_W(1)/1e5; // km/s;
      a_Wout(2) = a_W(2)/1e5; // km/s;
      a_Wout(3) = a_W(3)/1e5; // km/s;
      a_Wout(4) = a_W(4)*1e12; // picodyne/cm*2
      a_Wout(5) = a_W(5)*1e6; // microGauss
      a_Wout(6) = a_W(6)*1e6; // microGauss
      a_Wout(7) = a_W(7)*1e6; // microGauss
    }, Wout[dit], W);
  }
  std::string output_file = ParseInputs::get_data_file_prefix() + "_" + std::to_string(iter);
  h5.writeMBLevel({"density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp"}, map, Wout, output_file);
  // h5.writeMBLevel({"density","Momr","Momt","Momp","E","Br","Bt","Bp"}, map, JUTemp, "U_sph_" + to_string(iter));
  // h5.writeMBLevel({"density","Momr","Momt","Momp","E","Br","Bt","Bp"}, map, a_JU, "JU_" + to_string(iter));
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

    if (dbx1.contains(index_cc))
    {
      for (int i = 0; i < NUMCOMPS; i++)
      {
        probed_values[i] = U_dim(index_cc, i);
        if (dbx0.contains(index_n0))
          probed_values[i] += (U_dim(index_n0, i) - U_dim(index_cc, i)) * (pt0 - (pt0_nearest + 0.5)) / (pt0_neigbor - (pt0_nearest + 0.5));
        if (dbx0.contains(index_n1))
          probed_values[i] += (U_dim(index_n1, i) - U_dim(index_cc, i)) * (pt1 - (pt1_nearest + 0.5)) / (pt1_neigbor - (pt1_nearest + 0.5));
        if (dbx0.contains(index_n2))
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

  outputFile.close();
}
