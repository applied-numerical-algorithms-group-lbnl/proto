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

size_t dim(size_t idx, sycl::id<DIM> orng, sycl::id<DIM> irng) {
	sycl::id<DIM> id;
	id[DIM-1] = idx%orng[DIM-1] + 1;
#if DIM>1
	id[DIM-2] = (idx/orng[DIM-1])%orng[DIM-2] + 1;
#if DIM>2
	id[0] = (idx/orng[2])/orng[1] + 1;
#endif
#endif
	size_t ret = id[DIM-1];
#if DIM>1
	ret += id[DIM-2]*irng[DIM-1];
#if DIM>2
	ret += id[0]*irng[DIM-1]*irng[DIM-2];
#endif
#endif
	return ret;
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

template<typename T, size_t C=1>
inline void setup(T* (&ptrs)[C], T *buff, size_t size) {
	ptrs[0] = buff;
	for (size_t c=1; c<C; c++)
		ptrs[c] = ptrs[c-1]+size;
}

const sycl::range<DIM> range(const Box &bx) { 
#if DIM==1
	return sycl::range<1>(bx.size(0));
#elif DIM==2
	return sycl::range<2>(bx.size(0),bx.size(1));
#elif DIM==3
	return sycl::range<3>(bx.size(0),bx.size(1),bx.size(2));
#endif
	}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreachIP(const Func &f, BufferData<T,C> &b, Args... arg) {
	T* ptrs[C]; setup(ptrs, b.buff, b.bx.size());
	b.queue().parallel_for(b.range(),[=](sycl::item<DIM> item){
		Idx<T,C> idx(ptrs, item.get_linear_id());
		f(idx,arg...); 
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreachIP(const Func &f, BufferData<T,C> &in, BufferData<T,C> &out, Args... arg) {
	T* iptr[C]; setup(iptr, in.buff, in.bx.size());
	T* optr[C]; setup(optr, out.buff, out.bx.size());
	out.queue().parallel_for(range(in.bx&out.bx),[=](sycl::item<DIM> item){ 
		size_t id = item.get_linear_id();
		Idx<T,C> idx(iptr,id);
		Idx<T,C> odx(optr,id);
		f(idx,odx,arg...); 
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach(const Func &f, BufferData<T,C> &in, Args... arg) {
	BufferData<T,C> out(in.box());
	T* iptr[C]; setup(iptr, in.buff, in.bx.size());
	T* optr[C]; setup(optr, out.buff, out.bx.size());
	out.queue().parallel_for(out.range(),[=](sycl::item<DIM> item){ 
		size_t id = item.get_linear_id();
		Idx<T,C> idx(iptr,id);
		Idx<T,C> odx(optr,id);
		f(idx,odx,arg...); 
	});
	return out;
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach(const Func &f, BufferData<T,C> &left, BufferData<T,C> &right, Args... arg) {
	BufferData<T,C> ret(left.box()&right.box());
	T* iptr[C]; setup(iptr, left.buff, left.bx.size());
	T* optr[C]; setup(optr, right.buff, right.bx.size());
	T* rptr[C]; setup(rptr, ret.buff, ret.bx.size());
	ret.queue().parallel_for(ret.range(),[=](sycl::item<DIM> item){ 
		size_t id = item.get_linear_id();
		Idx<T,C> idx(iptr,id);
		Idx<T,C> odx(optr,id);
		Idx<T,C> rdx(rptr,id);
		f(idx,odx,rdx,arg...); 
	});
	return ret;
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreach_pIP(const Func &f, BufferData<T,C> &b, Args... arg) {
	sycl::id<DIM> base = util(b.box().low());
	T* ptr[C]; setup(ptr, b.buff, b.bx.size());
	b.queue().parallel_for(b.range(),[=](sycl::item<DIM> item){ 
		auto pt = base+item.get_id(); 
		Idx<T,C> idx(ptr,item.get_linear_id());
		f(idx,pt,arg...);
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline void foreach_pIP(const Func &f, BufferData<T,C> &in, BufferData<T,C> &out, Args... arg) {
	Box ixn = in.bx&out.bx;
	sycl::id<DIM> base = util(ixn.low());
	T* iptr[C]; setup(iptr, in.buff, in.bx.size());
	T* optr[C]; setup(optr, out.buff, out.bx.size());
	out.queue().parallel_for(range(ixn),[=](sycl::item<DIM> item){ 
		auto pt = base+item.get_id(); 
		size_t id = item.get_linear_id();
		Idx<T,C> idx(iptr,id);
		Idx<T,C> odx(optr,id);
		f(idx,odx,pt,arg...); 
	});
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach_p(const Func &f, BufferData<T,C> &b, Args... arg) {
	BufferData<T,C> ret(b.box());
	sycl::id<DIM> base = util(ret.bx.low());
	T* iptr[C]; setup(iptr, b.buff, b.bx.size());
	T* optr[C]; setup(optr, ret.buff, ret.bx.size());
	ret.queue().parallel_for(ret.range(),[=](sycl::item<DIM> item){ 
		auto pt = base+item.get_id(); 
		size_t id = item.get_linear_id();
		Idx<T,C> idx(iptr,id);
		Idx<T,C> odx(optr,id);
		f(idx,odx,pt,arg...); 
	});
	return ret;
}

template<typename Func, typename T, size_t C, typename... Args>
inline BufferData<T,C> foreach_p(const Func &f, BufferData<T,C> &left, BufferData<T,C> &right, Args... arg) {
	BufferData<T,C> ret(left.box()&right.box());
	sycl::id<DIM> base = util(ret.box().low());
	T* iptr[C]; setup(iptr, left.buff, left.bx.size());
	T* optr[C]; setup(optr, right.buff, right.bx.size());
	T* rptr[C]; setup(rptr, ret.buff, ret.bx.size());
	ret.queue().parallel_for(ret.range(),[=](sycl::item<DIM> item){ 
		auto pt = base+item.get_id(); 
		size_t id = item.get_linear_id();
		Idx<T,C> idx(iptr,id);
		Idx<T,C> odx(optr,id);
		Idx<T,C> rdx(rptr,id);
		f(idx,odx,rdx,pt,arg...);
	});
	return ret;
}

template <typename T, size_t C>
struct BufferData {
	friend class SimpleStencil;
	sycl::queue q;
	T* buff;
	Box bx;
	size_t size;
	BufferData() = delete;
	BufferData(const Box &box) : q{}, bx(box), size(C*box.size()), buff(sycl::malloc_device<T>(C*box.size(),q)) {
		assert(q.get_device().is_gpu()); } 

	BufferData(const Box &box, T val) : BufferData(box) {
		T* ptr = buff;
		q.parallel_for(size,[=](sycl::id<1> id){
			ptr[id] = val;
		}).wait();
	}
	
	sycl::queue& queue() { return q; }

	const sycl::range<DIM> range() const { 
#if DIM==1
		return sycl::range<1>(bx.size(0));
#elif DIM==2
		return sycl::range<2>(bx.size(0),bx.size(1));
#elif DIM==3
		return sycl::range<3>(bx.size(0),bx.size(1),bx.size(2));
#endif
	}

	const Box box() const { return bx; }

	void print() { 
		sycl::queue hq(sycl::cpu_selector_v);
		T* ptr = sycl::malloc_host<T>(size,hq);
		q.copy<T>(buff,ptr,size).wait();
		sycl::range<DIM> r = range();
		for (size_t c=0; c<C; c++) {
			sycl::buffer<T,DIM> pbuf(ptr,r);
			sycl::host_accessor<T,DIM,sycl::access_mode::read> acc(pbuf);
			for (int i=0; i<box().size(0); i++) {
#if DIM>1
				for (int j=0; j<box().size(1); j++) {
#endif
#if DIM>2
					for (int k=0; k<box().size(2); k++) {
						sycl::id<3> id(i,j,k);
#elif DIM==2
						sycl::id<2> id(i,j);
#elif DIM==1
						sycl::id<1> id(i);
#endif
						cout << acc[id] << ' ';
					}
#if DIM>1
					cout << endl;
				}
#endif
#if DIM>2
				cout << endl;
			}
#endif
			ptr+=bx.size();
			cout << endl << "===============" << endl;
		}
	}

	void random(T lo, T hi) {
		static std::mt19937 mt_engine(std::random_device{}());
		static std::uniform_real_distribution<T> idist(lo,hi);
		sycl::queue hq(sycl::cpu_selector_v);
		T* ptr = sycl::malloc_host<T>(size,hq);
		for (size_t i=0; i<size; i++) ptr[i] = idist(mt_engine);
		q.copy<T>(ptr,buff,size).wait();
		free(ptr,hq);
	}
};

template<typename T, size_t C>
class Apply {
public:
	Apply(T* _in, T* _out, pair<int,T> *_sten, sycl::id<DIM> _ir, sycl::id<DIM> _or)
		: in(_in), out(_out), sten(_sten), irng(_ir), orng(_or) {}
	void operator() (sycl::id<1> id) const {
		size_t idx = dim(id,orng,irng);
		T val = in[idx];
		for (int i=0; i<2*DIM+1; i++) {
			int offset = sten[i].first;
			T coef = sten[i].second;
			val += in[idx+offset]*coef;
		}
		out[id] = val;
	}
private:
	T *in, *out;
	pair<int,T> *sten;
	sycl::id<DIM> irng, orng;
};

struct SimpleStencil {
	template<typename T, size_t C>
	static BufferData<T,C> apply(BufferData<T,C> &in) {
		Point lo(in.box().low()+Point::Ones());
		Point hi(in.box().high()-Point::Ones());
		Box bx(lo,hi);
		BufferData<T,C> out(bx);
		size_t length = in.range()[0]; // assumes symmetric
		pair<int,T> *coefs = sycl::malloc_device<pair<int,T>>(2*DIM+1,in.q);
		in.q.single_task([=](){
			coefs[0] = make_pair(0,-2*DIM);
			for (size_t i=1; i<2*DIM+1; i+=2) {
				coefs[i]  = make_pair(-pow(length,i/2),1);
				coefs[i+1]= make_pair( pow(length,i/2),1);
			}
		}).wait();
		size_t inset = in.box().size(), outset = out.box().size();
		sycl::id<DIM> irng = util(in.box().sizes()), orng = util(out.box().sizes());
		for (int c=0; c<C; c++) 
			in.q.parallel_for(outset,Apply<T,C>(in.buff+c*inset,out.buff+c*outset,coefs,irng,orng)).wait();
		return out;
	}
};
