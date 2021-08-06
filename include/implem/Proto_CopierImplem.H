//===============================================
// COPIER PUBLIC API

template<class OP>
Copier<OP>::Copier(Copier<OP>& a_rhs)
{
}

template<class OP>
Copier<OP>::~Copier()
{
}

template<class OP>
void
Copier<OP>::clear()
{
}

template<class OP>
Copier<OP>&
Copier<OP>::operator=(const Copier<OP>& a_rhs)
{
}

template<class OP>
bool
Copier<OP>::operator==(const Copier<OP>& a_rhs) const
{
}

template<class OP>
void
Copier<OP>::reverse()
{
}

template<class OP>
void
Copier<OP>::execute()
{
}

//===============================================
// COPIER PRIVATE FUNCTIONS

template<class OP>
void
Copier<OP>::makeItSo(){}

template<class OP>
void
Copier<OP>::makeItSoBegin(){}

template<class OP>
void
Copier<OP>::makeItSoLocal(){}

template<class OP>
void
Copier<OP>::makeItSoEnd(){}

template<class OP>
void
Copier<OP>::completePendingSends() const{}

template<class OP>
void
Copier<OP>::allocateBuffers(){}

template<class OP>
void
Copier<OP>::writeToSendBuffers(){}

template<class OP>
void
Copier<OP>::postSends() const{}

template<class OP>
void
Copier<OP>::readFromRecvBuffers(){}

template<class OP>
void
Copier<OP>::postRecvs() const{}

//===============================================
// MOTIONITEM, BUFFERENTRY, and COPYITERATOR FUNCTIONS

template<class OP>
bool
MotionItem::operator==(const MotionItem a_rhs) const
{

}

template<class OP>
bool
BufferEntry::operator<(const BufferEntry& a_rhs) const
{
}


template<class OP>
CopyIterator::CopyIterator(const Copier<OP>& a_copier, MotionType)
{
}

template<class OP>
const MotionItem&
CopyIterator::operator()() const
{}

template<class OP>
const MotionItem&
CopyIterator::operator[](size_t a_index) const{}

template<class OP>
void
CopyIterator::operator++(){}

template<class OP>
void
CopyIterator::ok(){}

template<class OP>
void
CopyIterator::reset(){}

template<class OP>
size_t
CopyIterator::size(){}

