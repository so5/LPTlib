/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef DSLIB_COMM_DATA_BLOCK_H
#define DSLIB_COMM_DATA_BLOCK_H

#include <mpi.h>
#include "LPT_LogOutput.h"
#include "PMlibWrapper.h"

namespace DSlib
{
//! DataBlockのDataとTime以外のメンバをまとめた構造体
struct CommDataBlockHeader
{
  long      BlockID;
  int       SubDomainID;
  REAL_TYPE Origin[3];
  int       OriginCell[3];
  int       BlockSize[3];
  REAL_TYPE Pitch[3];
  bool      has_cut;
};

//! @brief データブロックのヘッダ部をまとめた構造体とデータ領域ヘのポインタ、それぞれの転送用MPI_Request変数をまとめて保持するクラス
//!
//! サイズ指定のコンストラクタを呼ぶと、データ領域用のメモリ(Buff)を指定されたサイズで確保する
//! 受信側はBuffは、DSlib::AddCache()内でキャッシュにポインタを移動されBuffにはNULLが代入されるので
//! デストラクタ内のdeleteでは解放されない。
//! 逆に送信側ではこのdleteで全て解放される
class CommDataBlockManager
{
private:
  CommDataBlockManager();
  CommDataBlockManager(const CommDataBlockManager& obj);
  CommDataBlockManager& operator=(const CommDataBlockManager& obj);

public:
  CommDataBlockManager(int size)
  {
    const int vlen = 3;
    try
    {
      Buff   = new REAL_TYPE[size * vlen];
      d_cut  = new long long[size];
      Header = new CommDataBlockHeader;
    }
    catch (std::bad_alloc)
    {
      std::cerr << "couldn't allocate memory for communication buffer" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    //交点フラグにゴミが残っているとマズイので0クリアする
    for (size_t i = 0; i < size; i++)
    {
      d_cut[i] = 0;
    }
  }

  ~CommDataBlockManager()
  {
    delete Header;
    delete[] d_cut;
    d_cut = NULL;
    delete[] Buff;
    Buff = NULL;
  }

  bool Wait()
  {
    LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
    PM.start("LPT_MPI_Wait");
    MPI_Status status;
    if (MPI_SUCCESS != MPI_Wait(&Request0, &status))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("MPI_Wait for Data Failed");
    }
    if (MPI_SUCCESS != MPI_Wait(&Request1, &status))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("MPI_Wait for Header Failed");
    }
    if (MPI_SUCCESS != MPI_Wait(&Request2, &status))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("MPI_Wait for cut Failed");
    }
    PM.stop("LPT_MPI_Wait");
    return true;
  }

  bool Test()
  {
    LPT::PMlibWrapper& PM = LPT::PMlibWrapper::GetInstance();
    PM.start("LPT_MPI_Wait");
    MPI_Status status;
    int        flag0 = 0;
    int        flag1 = 0;
    int        flag2 = 0;
    if (MPI_SUCCESS != MPI_Test(&Request0, &flag0, &status))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("MPI_Test for Data Failed");
    }
    if (MPI_SUCCESS != MPI_Test(&Request1, &flag1, &status))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("MPI_Test for Header Failed");
    }
    if (MPI_SUCCESS != MPI_Test(&Request2, &flag2, &status))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("MPI_Test for cut Failed");
    }
    PM.stop("LPT_MPI_Wait");

    // MPI_Test後に flag == 0の時、転送は未完了
    return flag0 != 0 && flag1 != 0 && flag2 != 0;
  }

  REAL_TYPE*           Buff;     //!< データ送受信用のバッファ
  long long*           d_cut;    //!< cut配列送受信用のバッファ
  CommDataBlockHeader* Header;   //!< DataBlockのヘッダ部分送信用構造体
  MPI_Request          Request0; //!< データ配列の送受信につかうMPI_Request変数
  MPI_Request          Request1; //!< ヘッダ部の送受信につかうMPI_Request変数
  MPI_Request          Request2; //!< cut配列の送受信につかうMPI_Request変数
};
} // namespace DSlib
#endif
