#include "BufferData.H"

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
void point_bothF(Idx<double,2> in, Idx<double,2> &out, sycl::id<DIM> p, double val) { 
	for (int i=0; i<2; i++)
		out(i) = in(i)+p[0]-val; 
}
SYCL_KERNEL_END(point_bothF, point_both)

int main() {
	typedef BufferData<double,2> BDtype;
	Box bx = Box::Cube(5);
cout << "here" << endl;
	BDtype start(bx);
	start.random(0.,10.);
	cout << "start:" << endl;
	start.print(); 
	foreachIP(first,start,2.3);
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
/*	foreachIP(both,*bdvec[0],*bdvec[1],4.1);
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
