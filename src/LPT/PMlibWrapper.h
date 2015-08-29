/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef PMLIB_WRAPPER_H
#define PMLIB_WRAPPER_H
#ifdef USE_PMLIB
#include <cstdlib>
#include "PerfMonitor.h"
#endif
#include "MPI_Manager.h"
namespace LPT
{
//! @brief PMlibのインスタンスへの参照を一つ保持し、GetInstanceを経由したグローバルなアクセスを提供する
//
//  コンパイル時にUSE_PMLIBが指定されていなければ、GetInstance以外の全メソッドは何もしない
class PMlibWrapper
{
  //Singletonパターン
  /// Initialize()が呼ばれたかどうかのフラグ
  bool initialized;

private:
  PMlibWrapper() : PM(NULL)
  {
    initialized = false;
  }

  PMlibWrapper(const PMlibWrapper& obj);
  PMlibWrapper& operator=(const PMlibWrapper& obj);
  ~PMlibWrapper() {}
public:
  static PMlibWrapper& GetInstance()
  {
    static PMlibWrapper instance;
    return instance;
  }

  void Initialize(pm_lib::PerfMonitor* argPM)
  {
    PM = argPM;
#ifdef USE_PMLIB
    if (initialized) return;
    initialized = true;
    //排他区間
    PM->setProperties("LPT_Initialize", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_PrepareCalc", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_DestroyStartPoints", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_EmitParticle", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_DestroyParticle", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_MakeRequestQ", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_SortParticle", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_CalcNumComm", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_CommNumComm", pm_lib::PerfMonitor::COMM);
    PM->setProperties("LPT_PrepareComm", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_AlltoAllRequest", pm_lib::PerfMonitor::COMM);
    PM->setProperties("LPT_P2PRequest", pm_lib::PerfMonitor::COMM);
    PM->setProperties("LPT_CalcParticle", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_Discard_Cache", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_ExchangeParticleContainers", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_Post", pm_lib::PerfMonitor::CALC);
    PM->setProperties("LPT_FileOutput", pm_lib::PerfMonitor::CALC);

    // CalcParticleセクション内の詳細区間
    PM->setProperties("LPT_CommDataF2P", pm_lib::PerfMonitor::COMM, false);
    PM->setProperties("LPT_MPI_Test", pm_lib::PerfMonitor::COMM, false);
    PM->setProperties("LPT_MPI_Wait", pm_lib::PerfMonitor::COMM, false);
    PM->setProperties("LPT_AddCache", pm_lib::PerfMonitor::CALC, false);
    PM->setProperties("LPT_PP_Transport", pm_lib::PerfMonitor::CALC, false);
    PM->setProperties("LPT_MoveParticle", pm_lib::PerfMonitor::CALC, false);
    PM->setProperties("LPT_DelSendBuff", pm_lib::PerfMonitor::CALC, false);

//! memo
//PP_Transport::Calc()内で呼びだされるルーチンにも区間を設定したことがあったが
//1call あたり1e-6以下の実行時間となっておりタイマの解像度以下だったため測定には使えかった。
#endif
  }

private:
///PMlibのインスタンス
#ifdef USE_PMLIB
  pm_lib::PerfMonitor* PM;
#endif

public:
  //Public Method
  void start(const std::string& key)
  {
#ifdef USE_PMLIB
    PM->start(key);
#endif
  }

  void stop(const std::string& key, double flopPerTask = 0.0, unsigned iterationCount = 1)
  {
#ifdef USE_PMLIB
    PM->stop(key, flopPerTask, iterationCount);
#endif
  }
};
} //namespace LPT
#endif
