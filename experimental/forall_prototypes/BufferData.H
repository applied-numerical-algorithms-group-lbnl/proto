#include "Proto_Box.H"
#include <random>
#include <sycl/sycl.hpp>
using namespace std;
using namespace Proto;

template<typename T>
class SimpleStencil;

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
	void print() { sycl::host_accessor acc(buff); for (auto it:acc) cout << it << ' '; cout << endl; }
	void iota() {
		sycl::host_accessor<T,DIM> acc(buff);
		std::iota(acc.begin(),acc.end(),0);
	}
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
class SimpleStencil {
public:
	SimpleStencil(const std::vector<std::pair<int,T>> &vec) : m_coefs(vec) {}

	BufferData<T> apply(BufferData<T> &in) {
		Point lo(in.box().low()+Point::Ones());
		Point hi(in.box().high()-Point::Ones());
		Box bx(lo,hi);
		BufferData<T> out(bx);
		sycl::queue q(sycl::gpu_selector_v);
		typedef pair<int,T> vec_t;
		size_t width=m_coefs.size();
		sycl::buffer<vec_t> stencil(m_coefs.data(),sycl::range<1>(width));
		q.submit([&](sycl::handler &cgh) { 
			sycl::accessor<T,DIM,sycl::access_mode::read> in_acc(in.buff, cgh);
			sycl::accessor<T,DIM,sycl::access_mode::write> out_acc(out.buff, cgh);
			sycl::accessor<vec_t,1,sycl::access_mode::read> vec_acc(stencil, cgh);
			cgh.parallel_for(out.buff.get_range(),
				Apply<T,decltype(in_acc),decltype(out_acc),decltype(vec_acc)>(in_acc,out_acc,vec_acc));
		});
		return out;
	}
private:
	std::vector<std::pair<int,T>> m_coefs;
};
