/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef PPLIB_INTEGRATOR_H
#define PPLIB_INTEGRATOR_H

#include <iostream>
#include <cmath>

//forward declaration
namespace DSlib
{
class DataBlock;
}

namespace PPlib
{
//! @brief 4次のルンゲ=クッタ法による速度場の積分を行なう。
class Integrator
{
private:
  Integrator(const Integrator& obj);
  Integrator& operator=(const Integrator& obj);

public:
  Integrator(Interpolator* ip) : ptrIP(ip) {}
  //! @brief 4次ルンゲ=クッタ法による速度場の積分を行なう
  //! @param DataBlock  [in]   計算対象の粒子が存在するデータブロック
  //! @param t_step  [in]    ルンゲ=クッタ積分の時間刻み
  //! @param x_i     [inout] 粒子座標
  // int RKG(const DSlib::DataBlock& DataBlock, const double t_step, REAL_TYPE x_i[3]);
  int RKG(const DSlib::DataBlock& DataBlock, const double t_step, REAL_TYPE x_i[3], const std::string& id);

private:
  //! @brief ルンゲ=クッタ用の係数を計算する。
  //! @param DataBlock   [in]   計算対象の粒子が存在するデータブロック
  //! @param x_i         [in]   粒子座標
  //! @param func        [out]  ルンゲ=クッタ積分用の係数を計算した結果を格納する
  //  bool GetIntegrand(const DSlib::DataBlock& DataBlock, const REAL_TYPE x_i[3], REAL_TYPE func[3]);
  bool GetIntegrand(const DSlib::DataBlock& DataBlock, const REAL_TYPE x_i[3], REAL_TYPE func[3], const std::string& id);

  Interpolator* ptrIP; //!< 流速の補間に使うInterpolatorオブジェクトへのポインタ
};
} // namespace PPlib
#endif
