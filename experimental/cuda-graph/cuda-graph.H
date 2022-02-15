#pragma once
#ifdef PROTO_CUDA

// Aim: apply different pairs of stencil/box on one patch
template<typename T> class FusedStencilGraph
{
  using Stencil = Proto::Stencil<T>;
  using Box = Proto::Box;

  public:

  // constructor
  FusedStencilGraph() 
  {
    graphCreated = false;
  }


  // send data into device memory 
  // m_stencils must be filled before
  void copyInfo();

  std::vector<Stencil*>& getStencils();

  // Apply the stencil kernel on several stencils and boxes on a same in patch and out patch
  // Be carefull that m_stencil and *a_box are filled with the same order
  template <unsigned int C, unsigned char D, unsigned char E>
  void cudaGraph( const Proto::BoxData<T,C,D,E> & a_src,
                             Proto::BoxData<T,C,D,E> & a_dst,
                             Proto::Box         * a_box,
                             bool                a_initToZero,
                             T                   a_scale) ;
 
  private:

#ifdef PROTO_CUDA
  // Used to call the close function on each stencil
  void close();
#endif
  // Should have the same size
  std::vector<Stencil*> m_stencils; ///< to store the stencil pointers 
  const T** m_coefs;                ///< to store m_stencils[i]->devCoefs()
  int* m_ssize;                     ///< to store sizes of m_coefs and m_offsets
  const Point** m_offsets;          ///< to store m_sencils[i]->devOffsets
  bool graphCreated;
  cudaGraph_t m_graph;
  cudaGraphExec_t m_instance;

};

#ifdef PROTO_CUDA
  template<typename T>
  void FusedStencilGraph<T>::close()
  {
    for(std::size_t it = 0 ; it < m_stencils.size() ; it++)
	if(!m_stencils[it]->closed()) m_stencils[it]->close(); // copy information to the device
  }
#endif

  template<typename T>
  void FusedStencilGraph<T>::copyInfo()
  {
    close();
  }
  // Return a reference on the member m_stencil - Used to fill this field
  template<typename T>
  std::vector<Proto::Stencil<T>*>& FusedStencilGraph<T>::getStencils() 
  { 
	  return m_stencils; 
  }

  template <typename T> 
  template <unsigned int C, unsigned char D, unsigned char E>
  void FusedStencilGraph<T>::cudaGraph( const Proto::BoxData<T,C,D,E> & a_src,
                            Proto::BoxData<T,C,D,E>  & a_dst,
                            Proto::Box         * a_box,
                            bool                a_initToZero,
                            T                   a_scale) 
  {
    PR_TIME("cudaApply(Fused)");
    const int nStencil = m_stencils.size();
    // set value according to size boxes
   
    assert(nStencil>0); 

    int blocks; 
    int strides; 


#ifdef PROTO_CUDA


    protoStream_t stream = Proto::DisjointBoxLayout::getCurrentStream();

    if(!graphCreated)
    {
      cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
      for(int it = 0 ; it < m_stencils.size() ; it++)
      {
        Proto::Stencil<T> * ptr = m_stencils[it];
        for (int ee = 0; ee < E; ee++)
        {
          for (int dd = 0; dd < D; dd++)
          {
            blocks=a_box[it].flatten(0).size();
            strides=a_box[it].size(0);
            
            assert(blocks > 0);
            assert(strides> 0);

            const T* srcPtr = a_src.dataPtr(0, dd, ee);
            T* dstPtr = a_dst.dataPtr(0, dd, ee);
            int k;
#if DIM==3
            k = a_box[it].size(2);
#else
            k = 1;
#endif
            // launch my kernel
            auto devCoefs  = ptr->devCoefs();
            auto offsets   = ptr->offsets().data();
            auto destShift = ptr->destShift();
            auto srcRatio  = ptr->srcRatio();
            auto destRatio = ptr->destRatio();
	    auto& coefs     = ptr->coefs();
//#define printDebug
#ifdef printDebug
 	    std::cout << " Iterator: " << it << std::endl;
 	    std::cout << " number of Coefs: " << coefs.size() << std::endl;
 	    std::cout << " a_scale: " << a_scale << std::endl;
 	    std::cout << " a_dst: box -> " << a_dst.box() << std::endl;
 	    std::cout << " a_src: box -> " << a_src.box() << std::endl;
 	    std::cout << " a_box: box -> " << a_box[it] << std::endl;
 	    std::cout << std::endl;
#endif
	    Proto::stencilIndexer3<T,C> <<<blocks, strides, 0, stream>>>(
			a_box[it].size(0),  a_box[it].size(1),
			dstPtr, srcPtr,
			a_src.box().size(), a_dst.box().size(),
			a_src.box().low(),  a_src.box().high(),
			a_dst.box().low(),  a_dst.box().high(),
			a_box[it].low(),  a_box[it].high(),
			a_initToZero,
			coefs.size(), devCoefs, offsets, k,
			destShift, srcRatio,  destRatio, a_scale);

          }
        }
      }

      cudaStreamEndCapture(stream, &m_graph);
      cudaGraphInstantiate(&m_instance, m_graph, NULL, NULL, 0);
      this->graphCreated=true;
    }
    cudaGraphLaunch(m_instance, stream);
    cudaStreamSynchronize(stream);

#else
    std::cout << " the CPU version of cudaApplyFused is not implemented " << std::endl;
    std::abort();
#endif 
    unsigned long long int c = 0;
    for(int it = 0 ; it < nStencil ; it++)
      c += m_stencils[it]->numFlops(a_box[it]);
	    
    PR_FLOPS(c);
#ifdef printDebug
    protoDeviceSynchronize();
    std::cout << " end apply graph " << std::endl;
#endif 
  }

#else

#endif
