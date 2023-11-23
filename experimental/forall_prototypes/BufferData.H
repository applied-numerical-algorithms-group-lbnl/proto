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

sycl::id<DIM> util(const Point pt) {
#if DIM==1
	sycl::id<1> base(pt[0]);
#elif DIM==2
	sycl::id<2> base(pt[0],pt[1]);
#elif DIM==3
	sycl::id<3> base(pt[0],pt[1],pt[2]);
#endif
	return base;
}

template<typename T, size_t C=1>
class BufferData;

template<typename T, size_t C=1>
struct Idx {
	Idx(T* const (&arr)[C], size_t offset) {
		for (size_t c=0; c<C; c++) buff[c] = arr[c]+offset;
	}

	T& operator()(size_t c=0) { return *(buff[c]); }

	T* buff[C];
};

//TODO: make one Idx, then offset into it when needed
template<typename Func, typename T, size_t C, typename... Args>
inline void foreachIP(const Func &f, BufferData<T,C> &b, Args... arg) {
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		T* ptrs[C];
		for (size_t c=0; c<C; c++) {
			auto acc = b.buff[c].template get_access<sycl::access_mode::read_write>(cgh);
			ptrs[c] = &(acc[0]);
		}
		cgh.parallel_for(b.range(),[=](sycl::item<DIM> it){
			Idx<T,C> idx(ptrs, it.get_linear_id());
			f(idx,arg...); 
		});
	});
}
/*
template<typename Func, typename T, size_t C, typename... Args>
inline void foreachIP(const Func &f, BufferData<T,C> &in, BufferData<T,C> &out, Args... arg) {
	Idx<T,C> idx{in.buffs()};
	Idx<T,C> odx{out.buffs()};
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto irng = idx.buff.get_range(), orng = odx.buf.get_range();
		auto r = min(irng,orng);
		cgh.parallel_for(r,[=](sycl::id<DIM> id){ f(idx,odx,arg...); });
	});

}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach(const Func &f, BufferData<T,C> &in, Args... arg) {
	BufferData<T,C> out(in.box());
	Idx<T,C> idx{in.buffs()};
	Idx<T,C> odx{out.buffs()};
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		cgh.parallel_for(idx.buff.get_range(),[=](sycl::id<DIM> id){ f(idx,odx,arg...); });
	});
	return out;
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach(const Func &f, BufferData<T,C> &left, BufferData<T,C> &right, Args... arg) {
	BufferData<T,C> ret(left.box()&right.box());
	Idx<T,C> idx{in.buffs()};
	Idx<T,C> odx{out.buffs()};
	Idx<T,C> rdx{ret.buffs()};
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		sycl::range
		cgh.parallel_for(r,[=](sycl::id<DIM> id){ f(lacc[id],racc[id],ret_acc[id],arg...); });
	});
	return ret;
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreach_pIP(const Func &f, BufferData<T,C> &b, Args... arg) {
	sycl::id<DIM> base = util(b.box().low());
	sycl::queue q(sycl::gpu_selector_v);
	q.submit([&](sycl::handler &cgh) {
		auto acc = b.get_buff().template get_access<sycl::access_mode::write>(cgh);
		cgh.parallel_for(acc.get_range(),[=](sycl::id<DIM> id){ 
			auto pt = base+id; f(acc[id],pt,arg...); });
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreach_pIP(const Func &f, BufferData<T,C> &in, BufferData<T,C> &out, Args... arg) {
	sycl::id<DIM> base = util(b.box().low());
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
	BufferData<T,C> ret(bx);
	sycl::id<DIM> base = util(bx.low());
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
	sycl::id<DIM> base = util(ret.box().low());
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
*/
template <typename T, size_t C>
struct BufferData {
	friend class SimpleStencil;
	vector<sycl::buffer<T,DIM>> buff;
	Box bx;
	BufferData() = delete;
#if DIM==1
	BufferData(const Box &box) : bx(box) { for (size_t c=0; c<C; c++) buff.push_back({sycl::range<1>(box.size(0))}); }
#elif DIM==2
	BufferData(const Box &box) : bx(box), buff(sycl::range<2>(box.size(0),box.size(1))) {}
#elif DIM==3
	BufferData(const Box &box) : bx(box), buff(sycl::range<3>(box.size(0),box.size(1),box.size(2))) {}
#endif

	BufferData(const Box &box, T val) : BufferData(box) {
		sycl::queue q(sycl::gpu_selector_v);
		q.submit([&](sycl::handler &cgh){
			for (int c=0; c<C; c++) {
				auto acc = buff[c].get_access<sycl::access_mode::write>(cgh);
				cgh.parallel_for(acc.get_range(),[=](sycl::id<DIM> id) {
					acc[id] = val;
				});
			}
		});
	}

	const size_t size() const { return C*bx.size(); }

	const sycl::range<DIM> range() const { return buff[0].get_range(); }

	const Box box() const { return bx; }

	void print(size_t c=0) { 
		assert(c<C);
		for (int c=0; c<C; c++) {
			sycl::host_accessor acc(buff[c]); 
			for (int i=0; i<bx.size(0); i++) {
#if DIM>1
				for (int j=0; j<bx.size(1); j++) {
#endif
#if DIM==3
					for (int k=0; k<bx.size(2); k++) {
						sycl::range r(i,j,k);
#elif DIM==2
						sycl::range r(i,j);
#elif DIM==1
						sycl::range r(i);
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
			cout << endl << "===============" << endl;
		}	
	}

	void random(T lo, T hi) {
		static std::mt19937 mt_engine(std::random_device{}());
		static std::uniform_real_distribution<T> idist(lo,hi);
		vector<T> *vec = new vector<T>(size());
		for (auto &el : *vec) el = idist(mt_engine);
		sycl::buffer<T,1> vbuf(vec->data(),vec->size());
		sycl::range<1> comp(size()/C);
		sycl::queue q(sycl::gpu_selector_v);
		for (int c=0; c<C; c++) {
			q.submit([&](sycl::handler &cgh){
				auto bacc = buff[c].template get_access<sycl::access_mode::write>(cgh);
				sycl::id<1> offset(c*comp.size());
				auto vacc = vbuf.template get_access<sycl::access_mode::read>(cgh,comp,offset);
				cgh.parallel_for(bacc.get_range(),[=](sycl::item<DIM> item){
					sycl::id<DIM> did = item.get_id();
					sycl::id<1> lid = item.get_linear_id();
					bacc[did] = vacc[lid];
				});
			});
		}
		delete vec;
	}

	decltype(buff[0]) get_buff(size_t c=0) { return buff[c]; }
};

template<typename T, size_t C, typename InAcc, typename OutAcc, typename Stencil>
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

struct SimpleStencil {
	template<typename T, size_t C>
	static BufferData<T,C> apply(BufferData<T,C> &in) {
		Point lo(in.box().low()+Point::Ones());
		Point hi(in.box().high()-Point::Ones());
		Box bx(lo,hi);
		BufferData<T,C> out(bx);
		size_t length = in.range()[0]; // assumes symmetric
		vector<pair<int,T>> coefs(2*DIM+1); // only using neighbors in each dimension
		coefs[0] = make_pair(0,-2*DIM);
		for (int i=1; i<coefs.size(); i+=2) {
			coefs[i]  = make_pair(-pow(length,i/2),1);
			coefs[i+1]= make_pair( pow(length,i/2),1);
		}
		typedef pair<int,T> vec_t;
		size_t width=coefs.size();
		sycl::buffer<vec_t> stencil(coefs.data(),sycl::range<1>(width));
		sycl::queue q(sycl::gpu_selector_v);
		for (int c=0; c<C; c++)
			q.submit([&](sycl::handler &cgh) { 
				auto in_acc = in.buff[c].template get_access<sycl::access_mode::read>(cgh);
				auto out_acc = out.buff[c].template get_access<sycl::access_mode::write>(cgh);
				auto vec_acc = stencil.template get_access<sycl::access_mode::read>(cgh);
				cgh.parallel_for(out_acc.get_range(),
					Apply<T,C,decltype(in_acc),decltype(out_acc),decltype(vec_acc)>(in_acc,out_acc,vec_acc));
			});
		return out;
	}
};
