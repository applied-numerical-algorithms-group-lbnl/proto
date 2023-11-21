#include "BufferData.H"

int main() {
	Box bx = Box::Cube(5);
	BufferData<double> start(bx);
	start.iota();
	start.print(); 
	foreach(first<double>(),start,2.3);
	start.print(); 
	foreach(other<double>(),start,2.3,1.2);
	start.print(); 
	SimpleStencil<double> easy;
	vector<BufferData<double>*> bdvec((bx.size(0)+1)/2);
	bdvec[0] = &start;
	for (int i=1; i<bdvec.size(); i++) 
		bdvec[i] = new BufferData<double>{easy.apply(*(bdvec[i-1]))}; 
	cout << "Initial buffer: " << endl;
	for (int i=0; i<bdvec.size(); i++) 
		bdvec[i]->print(); 
	foreach(both<double>(),*bdvec[0],*bdvec[1],4.1);
	bdvec[1]->print();
	foreach_p(both<double>(),*bdvec[1],*bdvec[0],3.2);
	bdvec[0]->print();
	return 0;
}
