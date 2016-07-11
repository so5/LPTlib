/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <iostream>
#include <vector>
#include <mpi.h>

#include "Communicator.h"
#include "DecompositionManager.h"
#include "DSlib.h"
#include "PMlibWrapper.h"
#include "LPT_LogOutput.h"

namespace DSlib
{
int Isend(double* buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request* request)
{
  return MPI_Isend(buf, count, MPI_DOUBLE, dest, tag, comm, request);
}

int Isend(float* buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request* request)
{
  return MPI_Isend(buf, count, MPI_FLOAT, dest, tag, comm, request);
}

int Irecv(double* buf, int count, int source, int tag, MPI_Comm comm, MPI_Request* request)
{
  return MPI_Irecv(buf, count, MPI_DOUBLE, source, tag, comm, request);
}

int Irecv(float* buf, int count, int source, int tag, MPI_Comm comm, MPI_Request* request)
{
  return MPI_Irecv(buf, count, MPI_FLOAT, source, tag, comm, request);
}

bool Communicator::CommRequest2(DSlib* ptrDSlib, std::list< CommDataBlockManager* >* RecvBuff, const int& fence)
{
  LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
  PM.start("LPT_P2PRequest");
  bool need_to_rerun   = false;
  int  Myrank          = LPT::MPI_Manager::GetInstance()->get_myrank_p();
  int  RecvBuffMemSize = 0;
  if (LPT::MPI_Manager::GetInstance()->is_particle_proc())
  {
    for (int rank_f = 0; rank_f < LPT::MPI_Manager::GetInstance()->get_nproc_f(); rank_f++)
    {
      std::vector< long >& queue       = *(ptrDSlib->RequestQueues.at(rank_f));
      int                  num_request = queue.size();
      if (num_request > MaxRequestSize)
      {
        //既定サイズを越えていたら戻り値をtrueに変える
        need_to_rerun = true;
        num_request   = MaxRequestSize;
      }
      if (num_request > 0)
      {
        int dst = LPT::MPI_Manager::GetInstance()->get_rank_f2w(rank_f);
        MPI_Put(&*(queue.begin()), num_request, MPI_LONG, dst, MaxRequestSize * Myrank, num_request, MPI_LONG, window);

        int tag = 0;
        for (int i = 0; i < num_request; ++i)
        {
          CommDataBlockManager* tmp = new CommDataBlockManager(MaxDataBlockSize);
          RecvBuffMemSize += MaxDataBlockSize;

          int ierr0 = Irecv(tmp->Buff, MaxDataBlockSize * 3, dst, tag++, MPI_COMM_WORLD, &(tmp->Request0));
          if (ierr0 != MPI_SUCCESS) LPT::LPT_LOG::GetInstance()->ERROR("Irecv for Data failed. Return value = ", ierr0);

          int ierr1 = MPI_Irecv(tmp->Header, 1, MPI_DataBlockHeader, dst, tag++, MPI_COMM_WORLD, &(tmp->Request1));
          if (ierr1 != MPI_SUCCESS) LPT::LPT_LOG::GetInstance()->ERROR("Irecv for Header failed. Return value = ", ierr1);

          int ierr2 = MPI_Irecv(tmp->d_cut, MaxDataBlockSize, MPI_LONG_LONG_INT, dst, tag++, MPI_COMM_WORLD, &(tmp->Request2));
          if (ierr2 != MPI_SUCCESS) LPT::LPT_LOG::GetInstance()->ERROR("Irecv for cut failedReturn value = ", ierr2);

          if (ierr0 == MPI_SUCCESS && ierr1 == MPI_SUCCESS && ierr2 == MPI_SUCCESS) RecvBuff->push_back(tmp);

          //リクエストを送信したブロックIDをDSlib::RequestedBlocksに登録
          ptrDSlib->AddRequestedBlocks(queue.at(i));
        }
        // 要求したブロックIDをDSlib::RequestQueuesから削除
        // num_requestの最大値はqueue.size()なので、第二引数がqueue.end()を越える可能性は無い
        queue.erase(queue.begin(), queue.begin() + num_request);
      }
    }
  }

  RecvBuffMemSize *= sizeof(REAL_TYPE);
  LPT::LPT_LOG::GetInstance()->LOG("Memory size for Recv Buffer = ", RecvBuffMemSize);

  LPT::LPT_LOG::GetInstance()->LOG("P2P request done");
  PM.stop("LPT_P2PRequest");
  return need_to_rerun;
}

void Communicator::SendDataBlock(REAL_TYPE* Data, int* Mask, const int& vlen, std::list< CommDataBlockManager* >* SendBuff)
{
  LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
  PM.start("LPT_CommDataF2P");

  int SendBuffMemSize = 0;

  //CommRequest2()で行なったデータブロックIDの通信を完了させる
  MPI_Win_fence(0, window);
  if (LPT::MPI_Manager::GetInstance()->is_fluid_proc())
  {
    for (int rank_p = 0; rank_p < LPT::MPI_Manager::GetInstance()->get_nproc_p(); rank_p++)
    {
      int   tag         = 0;
      long* BlockIDList = BlockIDsToSend + rank_p * MaxRequestSize;
      int   dst         = LPT::MPI_Manager::GetInstance()->get_rank_p2w(rank_p);

      //IDリストの先頭から-1が入っているところまでを順に読み取って、流速データのパッキング
      //MPI_Isendの送信、IDリストの初期化を行う
      for (int i = 0; i < MaxRequestSize; i++)
      {
        long BlockID = BlockIDList[i];
        if (BlockID == -1) break;
        BlockIDList[i]            = -1;
        CommDataBlockManager* tmp = new CommDataBlockManager(MaxDataBlockSize);
        SendBuffMemSize += MaxDataBlockSize;
        CommPacking(BlockID, Data, Mask, vlen, tmp->Buff, tmp->d_cut, tmp->Header);
        int SendSize = tmp->Header->BlockSize[0] * tmp->Header->BlockSize[1] * tmp->Header->BlockSize[2];
        int ierr0    = Isend(tmp->Buff, SendSize * 3, dst, tag++, MPI_COMM_WORLD, &(tmp->Request0));
        if (ierr0 != MPI_SUCCESS) LPT::LPT_LOG::GetInstance()->ERROR("Isend for Data failed Return value = ", ierr0);

        int ierr1 = MPI_Isend(tmp->Header, 1, MPI_DataBlockHeader, dst, tag++, MPI_COMM_WORLD, &(tmp->Request1));
        if (ierr1 != MPI_SUCCESS) LPT::LPT_LOG::GetInstance()->ERROR("Isend for Header failed Return value = ", ierr1);

        int ierr2 = MPI_Isend(tmp->d_cut, SendSize, MPI_LONG_LONG_INT, dst, tag++, MPI_COMM_WORLD, &(tmp->Request2));
        if (ierr2 != MPI_SUCCESS) LPT::LPT_LOG::GetInstance()->ERROR("Isend for cut failed Return value = ", ierr2);

        if (ierr0 == MPI_SUCCESS && ierr1 == MPI_SUCCESS && ierr2 == MPI_SUCCESS) SendBuff->push_back(tmp);
      }
    }
  }

  SendBuffMemSize *= sizeof(REAL_TYPE);
  LPT::LPT_LOG::GetInstance()->LOG("Memory size for Send Buffer = ", SendBuffMemSize);
  PM.stop("LPT_CommDataF2P");
}

void Communicator::CommPacking(const long& BlockID, REAL_TYPE* Data, int* Mask, const int& vlen, REAL_TYPE* SendBuff, long long* d_cut_buff, CommDataBlockHeader* Header)
{
  DecompositionManager* ptrDM  = DecompositionManager::GetInstance();
  int                   halo   = ptrDM->GetGuideCellSize();
  int                   MyRank = LPT::MPI_Manager::GetInstance()->get_myrank_f();
  if (!LPT::MPI_Manager::GetInstance()->is_fluid_proc())
  {
    //そもそもcomm_fluidに属していなければ呼ばれないはず
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  Header->BlockID       = BlockID;
  Header->SubDomainID   = MyRank;
  Header->BlockSize[0]  = ptrDM->GetBlockSizeX(BlockID) + 2 * halo;
  Header->BlockSize[1]  = ptrDM->GetBlockSizeY(BlockID) + 2 * halo;
  Header->BlockSize[2]  = ptrDM->GetBlockSizeZ(BlockID) + 2 * halo;
  Header->Pitch[0]      = ptrDM->Getdx();
  Header->Pitch[1]      = ptrDM->Getdy();
  Header->Pitch[2]      = ptrDM->Getdz();
  Header->Origin[0]     = ptrDM->GetBlockOriginX(BlockID);
  Header->Origin[1]     = ptrDM->GetBlockOriginY(BlockID);
  Header->Origin[2]     = ptrDM->GetBlockOriginZ(BlockID);
  Header->OriginCell[0] = ptrDM->GetBlockOriginCellX(BlockID);
  Header->OriginCell[1] = ptrDM->GetBlockOriginCellY(BlockID);
  Header->OriginCell[2] = ptrDM->GetBlockOriginCellZ(BlockID);

  int BlockLocalOffset = ptrDM->GetBlockLocalOffset(BlockID, MyRank);
  int SubDomainSize[3];
  SubDomainSize[0] = ptrDM->GetSubDomainSizeX(MyRank) + 2 * halo;
  SubDomainSize[1] = ptrDM->GetSubDomainSizeY(MyRank) + 2 * halo;
  SubDomainSize[2] = ptrDM->GetSubDomainSizeZ(MyRank) + 2 * halo;

  int indexS = 0;
  // 袖領域も含めて転送する
  for (int i = 0; i < vlen; i++)
  {
    for (int l = 0; l < Header->BlockSize[2]; l++)
    {
      for (int k = 0; k < Header->BlockSize[1]; k++)
      {
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (int j = 0; j < Header->BlockSize[0]; j++)
        {
          SendBuff[indexS++] = Data[BlockLocalOffset + DecompositionManager::Convert4Dto1D(j, k, l, i, SubDomainSize[0], SubDomainSize[1], SubDomainSize[2])] * Mask[BlockLocalOffset + DecompositionManager::Convert3Dto1D(j, k, l, SubDomainSize[0], SubDomainSize[1])];
        }
      }
    }
  }
  indexS = 0;
  for (int l = 0; l < Header->BlockSize[2]; l++)
  {
    for (int k = 0; k < Header->BlockSize[1]; k++)
    {
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
      for (int j = 0; j < Header->BlockSize[0]; j++)
      {
          d_cut_buff[indexS++] = d_cut[BlockLocalOffset + DecompositionManager::Convert3Dto1D(j, k, l, SubDomainSize[0], SubDomainSize[1])];
      }
    }
  }
}
} // namespace DSlib
