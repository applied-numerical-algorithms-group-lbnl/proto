#include "BufferData.H"

int main() {
	Box bx = Box::Cube(7);
	BufferData<double> start(bx);
	start.iota();
	vector<pair<int,double>> v(2*DIM+1,make_pair(0,0));
	v[0]=make_pair(0,-2*DIM);
	size_t edge = bx.size(0);
	for (int i=1; i<v.size(); i+=2) {
		v[i  ] = make_pair(-pow(edge,i/2),1);
		v[i+1] = make_pair(pow(edge,i/2),1);
	}
	for (auto it : v) cout << it.first << ',' << it.second << '|'; cout << endl;
	SimpleStencil<double> easy(v);
	vector<BufferData<double>*> bdvec(4);
	bdvec[0] = &start;
	for (int i=1; i<4; i++) 
		bdvec[i] = new BufferData<double>{easy.apply(*(bdvec[i-1]))}; 
	cout << "Initial buffer: ";
	for (int i=0; i<4; i++)
		bdvec[i]->print();
	return 0;
}
