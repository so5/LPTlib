/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef DSLIB_COMMUNICATOR_H
#define DSLIB_COMMUNICATOR_H

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <mpi.h>

#include "CommDataBlock.h"
#include "DataBlock.h"

namespace DSlib
{
//forward declaration
struct Cache;
class DSlib;

//! 粒子-流体間の通信を行なうクラス
class Communicator
{
  Communicator(const Communicator& obj);
  Communicator& operator=(const Communicator& obj);

public:
  // Constructor
  Communicator(const int& argMaxRequestSize, const int& argMaxDataBlockSize, long long* d_cut) : MaxRequestSize(argMaxRequestSize),
                                                                                                 MaxDataBlockSize(argMaxDataBlockSize), d_cut(d_cut)
  {
    int NumProcs = LPT::MPI_Manager::GetInstance()->get_nproc_p();
    //TODO MPI_Type_structを使ってヘッダ+DataBlock+cut配列という形式にする？
    // char型の配列としてsizeof(DataBlock)+sizeof(REAL_TYPE)*MaxDataBlockSize+sizeof(long long)*MaxDataBlockSizeの領域をバッファとして用意して
    // memcopyで代入するか、各領域の先頭ポインタを計算して配列として扱えば一応まとめて送れるはず。
    //
    // MPI_Type_contiguous(sizeof(DataBlock), MPI_BYTE, &MPI_TypeDataBlock);
    // if(MPI_SUCCESS != MPI_Type_commit(&MPI_TypeDataBlock))
    // {
    //     throw "MPI_Type_commit failed !";
    // }

    MPI_Type_contiguous(sizeof(CommDataBlockHeader), MPI_BYTE, &MPI_DataBlockHeader);
    if (MPI_SUCCESS != MPI_Type_commit(&MPI_DataBlockHeader))
    {
      throw "MPI_Type_commit failed !";
    }

    size_t BlockIDsToSendSize = LPT::MPI_Manager::GetInstance()->is_fluid_proc() ? NumProcs * MaxRequestSize : 0;
    if (LPT::MPI_Manager::GetInstance()->is_fluid_proc())
    {
      BlockIDsToSend = new long[BlockIDsToSendSize];
      for (int i = 0; i < BlockIDsToSendSize; i++)
      {
        BlockIDsToSend[i] = -1;
      }
    }
    MPI_Win_create(BlockIDsToSend, BlockIDsToSendSize * sizeof(long), sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
    MPI_Win_fence(0, window);
  }

  //Destructor
  ~Communicator()
  {
    if (LPT::MPI_Manager::GetInstance()->is_fluid_proc())
    {
      delete[] BlockIDsToSend;
    }
    MPI_Win_free(&window);
  }

  //! データブロック転送のリクエストを行いつつMPI_Irecvを発行する
  bool CommRequest2(DSlib* ptrDSlib, std::list< CommDataBlockManager* >* RecvBuff, const int& fence);

  //! @brief *Dataが示す領域に保持されているデータから、BlockIDに相当するブロックのデータを取り出して、SendBuffにパッキングする
  //! @param BlockID  [in]  必要な領域のブロックID
  //! @param Data     [in]  流体ソルバーからもらってきた物理量を格納しているデータ領域へのポインタ
  //! @param Mask     [in]  流体ソルバーからもらってきた物理量のマスク(物理量が存在しないセルは0他は1)
  //! @param vlen     [in]  Dataの領域に格納されている物理量のベクトル長
  //! @param SendBuff [out] 送信バッファ
  //! @param d_cut_buff [out] cut情報用の送信バッファ
  //! @param CommDataBlockHeader [out] ヘッダ情報の送信バッファ
  void CommPacking(const long& BlockID, REAL_TYPE* Data, int* Mask, const int& vlen, REAL_TYPE* SendBuff, long long* d_cut_buff, CommDataBlockHeader* Header);

  //! 要求されたデータブロックを送信しつつRequestIDの受付領域を初期化する
  void SendDataBlock(REAL_TYPE* Data, int* Mask, const int& vlen, std::list< CommDataBlockManager* >* SendBuff);

private:
  MPI_Datatype MPI_DataBlockHeader; //!< DataBlockのヘッダデータ送受信用データタイプ
  MPI_Datatype MPI_TypeDataBlock;   //!< DataBlockの送受信用データタイプ
  long*        BlockIDsToSend;      //!< 自Rankから各Rankへ転送するデータブロックのIDを保持する領域
  size_t       MaxRequestSize;      //!< 1プロセスから同時に受け付ける最大ブロックID数
  int          MaxDataBlockSize;    //!< 最も大きいデータブロックに含まれるセル数(袖領域も含む)
  MPI_Win      window;              //!< ブロックIDの転送領域用MPI_Win変数
  long long*   d_cut;               //!< cut情報を保持する配列
};
} // namespace DSlib
#endif
