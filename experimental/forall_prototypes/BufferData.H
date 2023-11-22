#include "Proto_Box.H"
#include <random>
#include <sycl/sycl.hpp>
using namespace std;
using namespace Proto;

struct SimpleStencil;

#define SYCL_KERNEL_START inline
#define SYCL_KERNEL_END(local_name, app_name)                  \
struct struct_##local_name {   \
    template <typename... T> \
    inline void operator()(T&... args) const { local_name(args...);} \
    const char* myname = #app_name; \
};        \
static struct_##local_name app_name;

template<typename T, size_t C=1>
class BufferData;

template<typename Func, typename T, size_t C, typename... Args>
inline void foreachIP(const Func &f, BufferData<T,C> &b, Args... arg) {
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto acc = b.get_buff().template get_access<sycl::access_mode::read_write>(cgh);
		cgh.parallel_for(acc.get_range(),[=](sycl::id<DIM> id){ f(acc[id],arg...); });
	});

}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreachIP(const Func &f, BufferData<T,C> &in, BufferData<T,C> &out, Args... arg) {
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto in_acc = in.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto out_acc = out.get_buff().template get_access<sycl::access_mode::read_write>(cgh);
		cgh.parallel_for(in_acc.get_range(),[=](sycl::id<DIM> id){ f(in_acc[id],out_acc[id],arg...); });
	});

}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach(const Func &f, BufferData<T,C> &in, Args... arg) {
	BufferData<T,C> out(in.box());
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto in_acc = in.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto out_acc = out.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(in_acc.get_range(),[=](sycl::id<DIM> id){ f(in_acc[id],out_acc[id],arg...); });
	});
	return out;
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach(const Func &f, BufferData<T,C> &left, BufferData<T,C> &right, Args... arg) {
	BufferData<T,C> ret(left.box()&right.box());
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto lacc = left.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto racc = right.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto ret_acc = right.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(ret_acc.get_range(),[=](sycl::id<DIM> id){ f(lacc[id],racc[id],ret_acc[id],arg...); });
	});
	return ret;
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreach_pIP(const Func &f, BufferData<T,C> &b, Args... arg) {
	Point low = b.box().low();
#if DIM==1
	sycl::id<1> base(low[0]);
#elif DIM==2
	sycl::id<2> base(low[0],low[1]);
#elif DIM==3
	sycl::id<3> base(low[0],low[1],low[2]);
#endif
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto acc = b.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(acc.get_range(),[=](sycl::id<DIM> id){ 
			auto pt = base+id; f(acc[id],pt,arg...); });
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreach_pIP(const Func &f, BufferData<T,C> &in, BufferData<T,C> &out, Args... arg) {
	Point low = in.box().low();
#if DIM==1
	sycl::id<1> base(low[0]);
#elif DIM==2
	sycl::id<2> base(low[0],low[1]);
#elif DIM==3
	sycl::id<3> base(low[0],low[1],low[2]);
#endif
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto in_acc = in.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto out_acc = out.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(out_acc.get_range(),[=](sycl::id<DIM> id){ 
			auto pt = base+id; f(in_acc[id],out_acc[id],pt,arg...); });
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach_p(const Func &f, BufferData<T,C> &b, Args... arg) {
	Box bx = b.box();
	Point low = bx.low();
	BufferData<T,C> ret(bx);
#if DIM==1
	sycl::id<1> base(low[0]);
#elif DIM==2
	sycl::id<2> base(low[0],low[1]);
#elif DIM==3
	sycl::id<3> base(low[0],low[1],low[2]);
#endif
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto in_acc = b.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto out_acc = ret.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(in_acc.get_range(),[=](sycl::id<DIM> id){ 
			auto pt = base+id; f(in_acc[id],out_acc[id],pt,arg...); });
	});
	return ret;
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach_p(const Func &f, BufferData<T,C> &left, BufferData<T,C> &right, Args... arg) {
	BufferData<T,C> ret(left.box()&right.box());
	Point low = ret.box().low();
#if DIM==1
	sycl::id<1> base(low[0]);
#elif DIM==2
	sycl::id<2> base(low[0],low[1]);
#elif DIM==3
	sycl::id<3> base(low[0],low[1],low[2]);
#endif
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto lacc = left.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto racc = right.get_buff().template get_access<sycl::access_mode::read>(cgh);
		auto ret_acc = right.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(ret_acc.get_range(),[=](sycl::id<DIM> id){ 
			auto pt = base+id; f(lacc[id],racc[id],ret_acc[id],pt,arg...); });
	});
	return ret;
}

template <typename T, size_t C>
class BufferData {
	friend class SimpleStencil;
	sycl::buffer<array<T,C>,DIM> buff;
	Box bx;
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
			cgh.parallel_for(r,[=](sycl::id<DIM> id) {
				for (int i=0; i<C; i++) acc[id][i] = val;});
		});
	}
	const size_t size() const { return C*bx.size(); }
	const Box box() const { return bx; }
	void print(size_t cc=0) { 
		assert(cc<C);
		sycl::host_accessor acc(buff); 
		for (int cc=0; cc<C; cc++) {
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
					 cout << acc[r][cc] << ' ';
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
		cout << "===============" << endl;
		}	
	}
	void random(T lo, T hi) {
		sycl::host_accessor<array<T,C>,DIM> acc(buff);
		static std::mt19937 mt_engine(std::random_device{}());
		static std::uniform_real_distribution<T> idist(lo,hi);
		for (auto& el : acc) 
			for (int i=0; i<C; i++)
			  el[i] = idist(mt_engine);
	}
	decltype(buff) get_buff() { return buff; }
};

template<typename T, size_t C, typename InAcc, typename OutAcc, typename Stencil>
class Apply {
public:
	Apply(InAcc _in, OutAcc _out, Stencil _sten) : in(_in), out(_out), sten(_sten) {}
	void operator() (sycl::item<DIM> it) const {
		sycl::id<DIM> in_id(it);
		for (int i=0; i<DIM; i++) in_id[i]++;
		for (int c=0; c<C; c++) {
			T val=in[in_id][c];
			for (int i=0; i<sten.size(); i++) {
				int offset = sten[i].first;
				T coef = sten[i].second;
				val += in[in_id+offset][c]*coef;
			}
			sycl::id<DIM> id = it.get_id();
			out[id][c]=val; 
		}
	}
private:
	InAcc in; OutAcc out; Stencil sten;
};

struct SimpleStencil {
	template<typename T, size_t C>
	static BufferData<T,C> apply(BufferData<T,C> &in) {
		Point lo(in.box().low()+Point::Ones());
		Point hi(in.box().high()-Point::Ones());
		Box bx(lo,hi);
		BufferData<T,C> out(bx);
		size_t length = in.buff.get_range()[0]; // assumes symmetric
		vector<pair<int,T>> coefs(2*DIM+1); // only using neighbors in each dimension
		coefs[0] = make_pair(0,-2*DIM);
		for (int i=1; i<coefs.size(); i+=2) {
			coefs[i]  = make_pair(-pow(length,i/2),1);
			coefs[i+1]= make_pair( pow(length,i/2),1);
		}
		sycl::queue q(sycl::gpu_selector_v);
		typedef pair<int,T> vec_t;
		size_t width=coefs.size();
		sycl::buffer<vec_t> stencil(coefs.data(),sycl::range<1>(width));
		q.submit([&](sycl::handler &cgh) { 
			auto in_acc = in.buff.template get_access<sycl::access_mode::read>(cgh);
			auto out_acc = out.buff.template get_access<sycl::access_mode::write>(cgh);
			auto vec_acc = stencil.template get_access<sycl::access_mode::read>(cgh);
			cgh.parallel_for(out.buff.get_range(),
				Apply<T,C,decltype(in_acc),decltype(out_acc),decltype(vec_acc)>(in_acc,out_acc,vec_acc));
		});
		return out;
	}
};
