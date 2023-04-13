
StencilTP<T>::build1D5pt(array<double,5>& a_coefs, int a_dir)
{
    return
          a_coefs[0]*Shift::Basis(a_dir,-2) +
          a_coefs[1]*Shift::Basis(a_dir,-1) +
          a_coefs[2]*Shift::Zeros() +
          a_coefs[3]*Shift::Basis(a_dir,1) +
          a_coefs[4]*Shift::Basis(a_dir,2);
}

// Construct tensor product stencils.
StencilTP<T>::Build(int a_order,int a_refratio)
{
    PROTO_ASSERT((a_order == 4),“Tensor Product Stencil Order must be 4”);
    array<array<T,5> , 4 > coefs4thRef4 =
    {{
         {-35.0/1024,294.0/1024,960.0/1024,-230.0/1024,35.0/1024},
         {-13.0/1024,58.0/1024,1088.0/1024,-122.0/1024,13.0/1024},
         {13.0/1024,-122.0/1024,1088.0/1024,58.0/1024,-13.0/1024},
         {35.0/1024,-230.0/1024,960.0/1024,294.0/1024,-35.0/1024}
     }};
    array<array<T,5> , 2 > coefs4thRef2 =
    {{
         {-3.0/128,22.0/128,128.0/128,-22.0/128,3.0/128},
         {3.0/128,-22.0/128,128.0/128,22.0/128,-3.0/128}
     }};
    array<InterpStencil<T>, DIM> m_cfInterpTP;
    for (int dir = 0; dir < DIM; dir++)
    {
        auto & cfidir = m_cfInterpTP[dir];
        cfidir.define(Point::Ones() + Point::Basis(dir,PR_AMR_REFRATIO-1));
        if (a_refratio == 4)
        {
            for (int k = 0; k < 4; k++)
                cfidir(Point::Basis(dir,k)) = build1D5pt(coefs4thRef4[k],dir);
        }
        else if (a_refratio == 2)
        {
            for (int k = 0; k < 2;k++)
                cfidir(Point::Basis(dir,k)) = build1D5pt(coefs4thRef2[k],dir);
        }
        else
        {
            PR_assert(false,“TP Stencil refratio must be 2 or 4");
        }
        cout << “size of stencil = ” << cfidir.size() << endl;
    }
}
StencilTP<T>::apply(BoxData<T,C,MEM>& a_input
        BoxData<T,C,MEM>& a_output)
{
    BoxData<T,C,MEM> UStageBox0(a_input.box());
    cout << “initial box = ” << a_input.box() << endl;
  a_input.copyTo(UStageBox0);
  for (int dir = 0; dir < DIM; dir++)
    {
      BoxData<T,C,MEM> UStageBox1 = m_cfInterpTP[dir](UStageBox0);
      cout << “dir = ” << dir << “, box = ” << UStageBox1.box() << endl;
      UStageBox0.define(UStageBox1.box());
      UStageBox1.copyTo(UStageBox0);
    }
  UStageBox0.copyTo(a_output);
}
