#include "BufferData.H"

int main() {
	Box bx = Box::Cube(5);
	BufferData<double> start(bx);
	start.iota();
	SimpleStencil<double> easy;
	vector<BufferData<double>*> bdvec((bx.size(0)+1)/2);
	bdvec[0] = &start;
	for (int i=1; i<bdvec.size(); i++) 
		bdvec[i] = new BufferData<double>{easy.apply(*(bdvec[i-1]))}; 
	cout << "Initial buffer: " << endl;
	for (int i=0; i<bdvec.size(); i++)
		bdvec[i]->print();
	return 0;
}
