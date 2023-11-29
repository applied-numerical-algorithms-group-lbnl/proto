#include "BufferData.H"

template<size_t C>
SYCL_KERNEL_START
void f_testF(const Point a_pt, const Box a_bx, Idx<double,C>& a_idx, const Array<double,C> &a_arr) {
	for (size_t c=0; c<C; c++) 
		a_idx(c) = a_bx.index(a_pt) * a_arr[C];
}
SYCL_KERNEL_END(f_testF, f_test)

template<size_t C>
   SYCL_KERNEL_START 
void f_eulerInitializeF(
        Point a_pt,
        Idx<double,C>& a_U,
        double a_h,
        double a_gamma,        
        double a_offset)
{
    double r0 = .125;
    double rsq = 0.;
    double amp = .01;
    for (int dir = 0; dir < DIM ; dir++)
    {
        double xcen = fmod(.5 + a_offset,1.); 
        double xdir = a_pt[dir]*a_h + .5*a_h;
        double del; 
        if (xcen > xdir) del = min(abs(xdir - xcen),abs(xdir - xcen + 1.));
        if (xcen <= xdir) del = min(abs(xdir - xcen),abs(xdir - xcen - 1.));
        rsq += pow(del,2);
    }
    double r = sqrt(rsq);
    if (r > r0)
    {
        a_U(0) = 1.;
    }
    else
    {
        a_U(0) = 1.0 + amp*pow(cos(M_PI*r/r0/2),6);
    }
    double ke = 0.;
    for (int dir = 0; dir < DIM; dir++)
    {
        a_U(dir+1) = 0.;
        ke += a_U(dir+1)*a_U(dir+1)/2;
        a_U(dir+1) *= a_U(0);
    }
    double pressure = pow(a_U(0),a_gamma);
    a_U(C-1) = a_gamma*pressure/(a_gamma - 1.) + a_U(0)*ke;
}
SYCL_KERNEL_END(f_eulerInitializeF, f_eulerInitialize)

SYCL_KERNEL_START
void firstF(Idx<double,2> &init, double val) {
	for (int i=0; i<2; i++)
		init(i) += val; 
}
SYCL_KERNEL_END(firstF, first)

SYCL_KERNEL_START
void otherF(Idx<double,2> &init, double first, double rest) { 
	for (int i=0; i<2; i++)
		init(i) = first+rest; 
}
SYCL_KERNEL_END(otherF, other)

SYCL_KERNEL_START
void bothF(Idx<double,2> in, Idx<double,2> &out, double val) { 
	for (int i=0; i<2; i++)
		out(i) = in(i)+val; 
}
SYCL_KERNEL_END(bothF, both)

SYCL_KERNEL_START
void allF(Idx<double,2> in, Idx<double,2> &out, double a, double b) { 
	for (int i=0; i<2; i++)
		out(i) = in(i)+a-b; 
}
SYCL_KERNEL_END(allF, all)

SYCL_KERNEL_START
void point_bothF(Idx<double,2> in, Idx<double,2> &out, const Point p, double val) { 
	for (int i=0; i<2; i++)
		out(i) = in(i)+p[0]-val; 
}
SYCL_KERNEL_END(point_bothF, point_both)

int main() {
	typedef BufferData<double,3> BDtype;
	Box bx = Box::Cube(5);
	BDtype start(bx);
	Array<double,3> arr{1.2,2.3,3.4};
	foreachTest(f_test,start,arr);
	//foreach_pIP(f_eulerInitialize,start,1.2,2.3,3.4);
//	start.random(0.,10.);
	cout << "start:" << endl;
	start.print(); 
/*	foreachIP(first,start,2.3);
	cout << "added 2.3:" << endl;
	start.print(); 
	foreachIP(other,start,2.3,1.2);
	cout << "set to 3.5:" << endl;
	start.print(); 
	vector<BDtype*> bdvec((bx.size(0)+1)/2);
	bdvec[0] = &start;
	for (int i=1; i<bdvec.size(); i++) 
		bdvec[i] = new BDtype{SimpleStencil::apply(*(bdvec[i-1]))}; 
	cout << "Initial buffer: " << endl;
	for (int i=0; i<bdvec.size(); i++) 
		bdvec[i]->print(); 
	foreachIP(both,*bdvec[0],*bdvec[1],4.1);
	cout << "adding 4.1:" << endl;
	bdvec[1]->print();
	cout << "adding vec[1]+p[0]-3.2 to vec[0]" << endl;
	foreach_pIP(point_both,*bdvec[1],*bdvec[0],3.2);
	bdvec[0]->print();
	auto last = foreach(all,start,2.4,1.3);
	cout << "final print:" << endl;
	last.print();
*/	return 0;
}
