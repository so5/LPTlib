/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <iostream>
#include <algorithm>

#include "DSlib.h"
#include "MPI_Manager.h"
#include "LPT_LogOutput.h"
#include "PMlibWrapper.h"

namespace DSlib
{
void DSlib::DiscardCacheEntry2(const long& num_entry)
{
  LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
  PM.start("LPT_PrepareComm");
  int room = CacheSize - CachedBlocks.size();
  if (room <= num_entry)
  {
    PurgeCachedBlocks(num_entry - room);
  }
  LPT::LPT_LOG::GetInstance()->LOG("DiscardCache done");
  PM.stop("LPT_PrepareComm");
}

void DSlib::AddRequestQueues(const int& SubDomainID, const long& BlockID)
{
  RequestQueues.at(SubDomainID)->push_back(BlockID);
}

long DSlib::AddCachedBlocks(CommDataBlockManager* RecvData, const double& Time)
{
  LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
  PM.start("LPT_AddCache");
  int ArrivedBlockID = -1;
  omp_set_lock(&CachedBlocksLock);
  if (CachedBlocks.size() > CacheSize)
  {
    LPT::LPT_LOG::GetInstance()->WARN("DataBlock Cache overflowed");
    LPT::LPT_LOG::GetInstance()->WARN("CachedBlocks size = ", CachedBlocks.size());
    LPT::LPT_LOG::GetInstance()->WARN("Max cache size = ", CacheSize);
    return -1;
  }
  omp_unset_lock(&CachedBlocksLock);

  ArrivedBlockID   = RecvData->Header->BlockID;
  DataBlock* tmp   = new DataBlock;
  tmp->BlockID     = ArrivedBlockID;
  tmp->SubDomainID = RecvData->Header->SubDomainID;
  tmp->Time        = Time;
  for (int i = 0; i < 3; i++)
  {
    tmp->Origin[i]     = RecvData->Header->Origin[i];
    tmp->OriginCell[i] = RecvData->Header->OriginCell[i];
    tmp->BlockSize[i]  = RecvData->Header->BlockSize[i];
    tmp->Pitch[i]      = RecvData->Header->Pitch[i];
  }
  tmp->Data       = RecvData->Buff;
  RecvData->Buff  = NULL;
  tmp->d_cut      = RecvData->d_cut;
  RecvData->d_cut = NULL;

  Cache* tmp2   = new Cache;
  tmp2->BlockID = ArrivedBlockID;
  tmp2->ptrData = tmp;
  for (size_t i = 0; i < tmp->BlockSize[0] * tmp->BlockSize[1] * tmp->BlockSize[2]; i++)
  {
    if ((tmp->d_cut)[i] & 0x3f != 0)
    {
      tmp->has_cut = true;
      break;
    }
  }

  omp_set_lock(&CachedBlocksLock);
  CachedBlocks.push_back(tmp2);
  omp_unset_lock(&CachedBlocksLock);

  PM.stop("LPT_AddCache");
  return ArrivedBlockID;
}

void DSlib::PurgeCachedBlocks(const int& NumEntry)
{
  if (NumEntry >= CachedBlocks.size())
  {
    omp_set_lock(&CachedBlocksLock);
    for (std::CACHE_CONTAINER< Cache* >::iterator it = CachedBlocks.begin(); it != CachedBlocks.end(); ++it)
    {
      delete (*it)->ptrData;
      delete *it;
    }
    std::CACHE_CONTAINER< Cache* >().swap(CachedBlocks);
    omp_unset_lock(&CachedBlocksLock);
    LPT::LPT_LOG::GetInstance()->LOG("All CachedBlocks is purged");
    return;
  }
  int DeleteCount = 0;

  std::CACHE_CONTAINER< Cache* >::iterator FirstEntry;
  for (std::CACHE_CONTAINER< Cache* >::reverse_iterator rit = CachedBlocks.rbegin(); rit != CachedBlocks.rend();)
  {
    if (DeleteCount == NumEntry)
    {
      break;
    }
    delete (*rit)->ptrData;
    delete *rit;
    FirstEntry = rit.base();
    DeleteCount++;
  }
  omp_set_lock(&CachedBlocksLock);
  CachedBlocks.erase(FirstEntry, CachedBlocks.end());
  omp_unset_lock(&CachedBlocksLock);
  std::CACHE_CONTAINER< Cache* >(CachedBlocks).swap(CachedBlocks);
}

void DSlib::PurgeAllCacheLists(void)
{
  LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
  PM.start("LPT_Discard_Cache");
  for (std::vector< std::vector< long >* >::iterator it = RequestQueues.begin(); it != RequestQueues.end(); ++it)
  {
    std::vector< long >().swap(**it);
  }
  RequestedBlocks.clear();
  PurgeCachedBlocks(CachedBlocks.size());
  PM.stop("LPT_Discard_Cache");
}

int DSlib::Load(const long& BlockID, DataBlock** DataBlock)
{
  //CachedBlocksの中を探索
  bool found = false;
  omp_set_lock(&CachedBlocksLock);
  for (std::CACHE_CONTAINER< Cache* >::iterator it = CachedBlocks.begin(); it != CachedBlocks.end(); ++it)
  {
    if ((*it)->BlockID == BlockID)
    {
      *DataBlock       = (*it)->ptrData;
      Cache* tmp_Cache = (*it);
      CachedBlocks.erase(it);
      CachedBlocks.push_front(tmp_Cache);
      found = true;
      break;
    }
  }
  omp_unset_lock(&CachedBlocksLock);
  if (found) return 0;

  //RequestedBlocksの中を探索
  if (RequestedBlocks.end() != RequestedBlocks.find(BlockID))
  {
    LPT::LPT_LOG::GetInstance()->LOG("DataBlock is not arrived: ", BlockID);
    return 1;
  }

  //RequestQueueの中を探索
  for (std::vector< std::vector< long >* >::iterator it = RequestQueues.begin(); it != RequestQueues.end(); ++it)
  {
    for (std::vector< long >::iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2)
    {
      if ((*it2) == BlockID)
      {
        LPT::LPT_LOG::GetInstance()->LOG("DataBlock is not requested at this time: ", BlockID);
        return 2;
      }
    }
  }
  //どこにも無い -> このタイムステップでの対象粒子の計算は中止
  LPT::LPT_LOG::GetInstance()->WARN("DataBlock is not requested at this time step: ", BlockID);
  return 4;
}
} // namespace DSlib
