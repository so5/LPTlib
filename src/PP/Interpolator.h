/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef PPLIB_INTERPOLATOR_H
#define PPLIB_INTERPOLATOR_H

//forward declaration
namespace DSlib
{
class DataBlock;
}
namespace PPlib
{
//! @brief 流速の補間を行なうクラス
class Interpolator
{
public:
  Interpolator(const Interpolator& arg)
  {
    for (int i = 0; i < 8; i++)
    {
      masks[i] = arg.masks[i];
    }
  }
  Interpolator& operator=(const Interpolator& arg)
  {
    for (int i = 0; i < 8; i++)
    {
      masks[i] = arg.masks[i];
    }
  }

  Interpolator()
  {
    masks[0] = 0x2aLL; //!< x plus,  y plus,  z plusの交点フラグを取り出すマスク
    masks[1] = 0x29LL; //!< x minus, y plus,  z plusの交点フラグを取り出すマスク
    masks[2] = 0x25LL; //!< x minus, y minus, z plusの交点フラグを取り出すマスク
    masks[3] = 0x26LL; //!< x plus,  y minus, z plusの交点フラグを取り出すマスク
    masks[4] = 0x1aLL; //!< x plus,  y plus,  z minusの交点フラグを取り出すマスク
    masks[5] = 0x19LL; //!< x minus, y plus,  z minusの交点フラグを取り出すマスク
    masks[6] = 0x15LL; //!< x minus, y minus, z minusの交点フラグを取り出すマスク
    masks[7] = 0x16LL; //!< x plus,  y minus, z minusの交点フラグを取り出すマスク
  }
  //! @brief ベクトルデータの補間を行う
  //! @param x_I  [in]  粒子座標
  //! @param didx [in]  補間対象データのindex番号
  //! @param dval [out] 補間したベクトルデータを格納する領域
  //  bool InterpolateData(const DSlib::DataBlock& DataBlock, const REAL_TYPE x_I[3], REAL_TYPE dval[3]);
  bool InterpolateData(const DSlib::DataBlock& DataBlock, const REAL_TYPE x_I[3], REAL_TYPE dval[3], const std::string id);

  //! @brief 交点を考慮したベクトルデータの補間を行う
  //! @param x_I  [in]  粒子座標
  //! @param didx [in]  補間対象データのindex番号
  //! @param dval [out] 補間したベクトルデータを格納する領域
  bool InterpolateWithCut(REAL_TYPE dval[3], REAL_TYPE data[8][3], const long long cuts[8], const REAL_TYPE& ip, const REAL_TYPE& jp, const REAL_TYPE& kp, const REAL_TYPE& im, const REAL_TYPE& jm, const REAL_TYPE& km);

  //! @brief 解析領域全体でのグローバル座標の座標値を、データブロック内のローカル座標に変換する
  //! @param [in]  x   解析領域内でのグローバル座標
  //! @param [out] x_i データブロック内でのローカル座標
  //! グローバル座標では、袖領域を含まない範囲で一番端のセルのコーナーを原点とする
  //! ローカル座標では、袖領域を含む範囲で一番端のセルの中心を原点とする
  void ConvXtoI(const REAL_TYPE x_g[3], REAL_TYPE x_l[3], const REAL_TYPE orig[3], const REAL_TYPE pitch[3]);

  //! @brief ConvXtoIの逆変換を行う
  //! @param x_I [in]  データブロック内でのローカル座標
  //! @param x   [out] 解析領域内でのグローバル座標
  void ConvItoX(const REAL_TYPE x_l[3], REAL_TYPE x_g[3], const REAL_TYPE orig[3], const REAL_TYPE pitch[3]);

private:
  //! @brief 与えられた8点の値を元にトリリニア補間を行う
  //! @param dval 補間後のベクトルデータを格納する領域
  //! @param data 補間に用いる周辺のデータ
  //! @param ip   補間点からx plus方向の既知の点までの距離
  //! @param jp   補間点からy plus方向の既知の点までの距離
  //! @param kp   補間点からz plus方向の既知の点までの距離
  //! @param im   補間点からx minus方向の既知の点までの距離
  //! @param jm   補間点からy minus方向の既知の点までの距離
  //! @param km   補間点からz minus方向の既知の点までの距離
  //
  bool TrilinearInterpolate(REAL_TYPE dval[3], REAL_TYPE data[8][3], const REAL_TYPE& ip, const REAL_TYPE& jp, const REAL_TYPE& kp, const REAL_TYPE& im, const REAL_TYPE& jm, const REAL_TYPE& km);

  long long masks[8]; //< cut情報から補間に関連する交点フラグを取り出すためのマスク
};
} // namespace PPlib
#endif
