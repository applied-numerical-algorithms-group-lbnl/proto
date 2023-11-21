#include "Proto_Box.H"
#include <random>
#include <sycl/sycl.hpp>
using namespace std;
using namespace Proto;

template<typename T>
class SimpleStencil;

template<typename T>
struct first {
	void operator()(T& init, T first) const { init += first; }
};

template<typename T>
struct other {
	void operator()(T& init, T first, T rest) const { init = first+rest; }
};

template<typename T>
struct both {
	void operator()(T in, T& out, T val) const { out = in+val; }
};

template<typename T>
class BufferData;

template<typename Func, typename T, typename... Args>
void foreach(const Func &f, BufferData<T> &b, Args... arg) {
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		sycl::accessor<T,DIM> acc = b.get_buff().get_access(cgh);
		cgh.parallel_for(b.get_buff().get_range(),[=](sycl::id<DIM> id){ f(acc[id],arg...); });
	});
}

template<typename Func, typename T, typename... Args>
void foreach(const Func &f, BufferData<T> &in, BufferData<T> &out, Args... arg) {
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		sycl::accessor<T,DIM> in_acc = in.get_buff().get_access(cgh);
		sycl::accessor<T,DIM> out_acc = out.get_buff().get_access(cgh);
		sycl::range<DIM> in_r = in_acc.get_range(), out_r = out_acc.get_range();
		sycl::range<DIM> r = (in_r.size() < out_r.size()) ? in_r : out_r;
		cgh.parallel_for(r,[=](sycl::id<DIM> id){ f(in_acc[id],out_acc[id],arg...); });
	});
}

template <typename T>
class BufferData {
	friend class SimpleStencil<T>;
public:
	BufferData() = delete;
#if DIM==1
	BufferData(const Box &box) : bx(box), buff(sycl::range<1>(box.size(0))) {}
#elif DIM==2
	BufferData(const Box &box) : bx(box), buff(sycl::range<2>(box.size(0),box.size(1))) {}
#elif DIM==3
	BufferData(const Box &box) : bx(box), buff(sycl::range<3>(box.size(0),box.size(1),box.size(2))) {}
#endif
	BufferData(const Box &box, T val) : BufferData(box) {
		sycl::range<DIM> r = buff.get_range();
		sycl::queue q(sycl::gpu_selector_v);
		q.submit([&](sycl::handler &cgh){
			auto acc = buff.get_access<sycl::access_mode::write>(cgh);
			cgh.parallel_for(r,[=](sycl::id<DIM> id){acc[id] = val;});
		});
	}
	const size_t size() const { return bx.size(); }
	const Box box() const { return bx; }
	void print() { 
		sycl::host_accessor acc(buff); 
#if DIM>2
		for (int i=0; i<bx.size(DIM-3); i++) {
#endif
#if DIM>1
			for (int j=0; j<bx.size(DIM-2); j++) {
#endif
#if DIM>0
				for (int k=0; k<bx.size(DIM-1); k++) {
#endif
#if DIM==3
					sycl::range r(i,j,k);
#elif DIM==2
					sycl::range r(j,k);
#elif DIM==1
					sycl::range r(k);
#endif
					 cout << acc[r] << ' ';
				}
#if DIM>1
				cout << endl;
			}
#endif
#if DIM>2
			cout << endl;
		}
#endif
		cout << endl;
	}
	void iota() {
		sycl::host_accessor<T,DIM> acc(buff);
		std::iota(acc.begin(),acc.end(),0);
	}
	sycl::buffer<T,DIM> get_buff() { return buff; }
private:
	sycl::buffer<T,DIM> buff;
	Box bx;
};

template<typename T, typename InAcc, typename OutAcc, typename Stencil>
class Apply {
public:
	Apply(InAcc _in, OutAcc _out, Stencil _sten) : in(_in), out(_out), sten(_sten) {}
	void operator() (sycl::item<DIM> it) const {
		sycl::id<DIM> in_id(it);
		for (int i=0; i<DIM; i++) in_id[i]++;
		T val=in[in_id];
		for (int i=0; i<sten.size(); i++) {
			int offset = sten[i].first;
			T coef = sten[i].second;
			val += in[in_id+offset]*coef;
		}
		sycl::id<DIM> id = it.get_id();
		out[id]=val; 
	}
private:
	InAcc in; OutAcc out; Stencil sten;
};

template<typename T>
struct SimpleStencil {
	BufferData<T> apply(BufferData<T> &in) {
		Point lo(in.box().low()+Point::Ones());
		Point hi(in.box().high()-Point::Ones());
		Box bx(lo,hi);
		BufferData<T> out(bx);
		size_t length = in.buff.get_range()[0]; // assumes symmetric
		vector<pair<int,T>> coefs(2*DIM+1); // only using neighbors in each dimension
		coefs[0] = make_pair(0,-2*DIM);
		for (int i=1; i<coefs.size(); i+=2) {
			coefs[i]  = make_pair(-pow(length,i/2),1);
			coefs[i+1]= make_pair( pow(length,i/2),1);
		}
		for (auto it : coefs) cout << it.first << ',' << it.second << '|'; cout << endl;
		sycl::queue q(sycl::gpu_selector_v);
		typedef pair<int,T> vec_t;
		size_t width=coefs.size();
		sycl::buffer<vec_t> stencil(coefs.data(),sycl::range<1>(width));
		q.submit([&](sycl::handler &cgh) { 
			sycl::accessor<T,DIM,sycl::access_mode::read> in_acc(in.buff, cgh);
			sycl::accessor<T,DIM,sycl::access_mode::write> out_acc(out.buff, cgh);
			sycl::accessor<vec_t,1,sycl::access_mode::read> vec_acc(stencil, cgh);
			cgh.parallel_for(out.buff.get_range(),
				Apply<T,decltype(in_acc),decltype(out_acc),decltype(vec_acc)>(in_acc,out_acc,vec_acc));
		});
		return out;
	}
};
