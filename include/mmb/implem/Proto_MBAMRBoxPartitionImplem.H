MBAMRBoxPartition::MBAMRBoxPartition(
        const MBProblemDomain&      a_coarseDomain,
        const std::vector<Point>&   a_refRatios)
{
    define(a_coarseDomain, a_refRatios);
}

void MBAMRBoxPartition::define(
        const MBProblemDomain&      a_coarseDomain,
        const std::vector<Point>&   a_refRatios)
{
    int maxLevels = a_refRatios.size()+1;
    MBProblemDomain domain = a_coarseDomain;
    for (int ii = 0; ii < maxLevels; ii++)
    {
        if (ii > 0)
        {
            domain = domain.refine(a_refRatios[ii-1]); 
        }
        m_levelPartitions.push_back(std::make_shared<MBBoxPartition>(domain));
    }
}

void MBAMRBoxPartition::loadBalance()
{
    for (int li = 0; li < numLevels(); li++)
    {
        // FIXME: this is temporary
        m_levelPartitions[li]->loadBalance();
    }
}

bool MBAMRBoxPartition::compatible(const MBAMRBoxPartition& a_rhs) const
{
    for (int li = 0; li < numLevels(); li++)
    {
        if (!m_levelPartitions[li]->compatible(*a_rhs.levelPartition(li))) {return false;}
    }
    return true;
}

std::shared_ptr<BoxPartition> MBAMRBoxPartition::blockPartition(
        unsigned int a_level,
        unsigned int a_block) const
{
    return m_levelPartitions[a_level]->blockPartition(a_block);
}
std::shared_ptr<MBBoxPartition> MBAMRBoxPartition::levelPartition(
        unsigned int a_level) const
{
    return m_levelPartitions[a_level];
}

void MBAMRBoxPartition::print()
{
    for (int ii = 0; ii < numLevels(); ii++)
    {
        m_levelPartitions[ii]->print();
    }
}
