/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <set>
#include <bitset>

#include "Interpolator.h"
#include "DecompositionManager.h"
#include "LPT_LogOutput.h"
#include "DataBlock.h"
#include "Utility.h"

namespace PPlib
{
  template <typename T>
  struct vector3
  {
    vector3(const T& arg_x, const T& arg_y, const T& arg_z) 
    {
      x=arg_x;
      y=arg_y;
      z=arg_z;
    }
    T x;
    T y;
    T z;

    bool operator < (const vector3& obj) const
    {
      if (this->x != obj.x) return this->x < obj.x;
      if (this->y != obj.y) return this->y < obj.y;
      return this->z < obj.z;
    }
  };

namespace
{
// 無名namespace内はFFVCのFB_Define.hから一部をコピーしたもの
#define MASK_9 0x1ff // 9 bit幅
#define QT_9 511     // 9bit幅の最大値
#define TOP_CUT 6

/*
     * @brief cut indexから指定方向の量子化値をとりだす
     * @param [in] c    cut index
     * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
     */
inline int getBit9(const long long c, const int dir)
{
  long long a = MASK_9;
  return (int)(((c >> TOP_CUT) >> dir * 9) & a);
}
/*
	 * @brief cut indexから指定方向交点の有無を返す
	 * @param [in] c    cut index
	 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
	 * @retval 交点あり(1)、なし(0)
	 */
inline int ensCut(const long long c, const int dir)
{
  long long a = 1;
  return (int)((c >> dir) & a);
}
/*
 * @brief 9bit幅の量子化
 * @param [in]  a  入力数値
 * @note -1/(2*511) < a < 1/(2*511)のとき、s=0
 *       a= 1.0 --> s=511
 *       0.0 <= a <= 1.0 を想定
 */

inline int quantize9(REAL_TYPE a)
{
  int       s;
  REAL_TYPE x = a * (REAL_TYPE)QT_9;

  if (x > 0.0)
  {
    s = (int)floor(x + 0.5);
  } else {
    s = (int)(-1.0 * floor(fabs(x) + 0.5));
  }

  if (s < 0 || QT_9 < s)
  {
    printf("quantize error in Interpolator : out of range %f > %d\n", a, s);
    exit(0);
  }

  return s;
}
#undef MASK_9
#undef QT_9
#undef TOP_CUT
}

bool Interpolator::TrilinearInterpolate(REAL_TYPE dval[3], REAL_TYPE data[8][3], const REAL_TYPE& ip, const REAL_TYPE& jp, const REAL_TYPE& kp, const REAL_TYPE& im, const REAL_TYPE& jm, const REAL_TYPE& km)
{
  for (int l = 0; l < 3; l++)
  {
    dval[l] = im * jm * km * data[0][l] + ip * jm * km * data[1][l] + ip * jp * km * data[2][l] + im * jp * km * data[3][l] + im * jm * kp * data[4][l] + ip * jm * kp * data[5][l] + ip * jp * kp * data[6][l] + im * jp * kp * data[7][l];
  }
  return true;
}

/*
   *
   *  cell number
   *
   *    7-----6
   *   /|    /|
   *  4-----5 |
   *  | 3---|-2
   *  |/    |/
   *  0-----1
   *
   */
bool Interpolator::InterpolateWithCut(REAL_TYPE dval[3], REAL_TYPE data[8][3], const long long cuts[8], const REAL_TYPE& ip, const REAL_TYPE& jp, const REAL_TYPE& kp, const REAL_TYPE& im, const REAL_TYPE& jm, const REAL_TYPE& km, const REAL_TYPE Origin[3], const REAL_TYPE Pitch[3])
{
  const int xdirs[8] = { 1, 0, 0, 1, 1, 0, 0, 1 };
  const int ydirs[8] = { 3, 3, 2, 2, 3, 3, 2, 2 };
  const int zdirs[8] = { 5, 5, 5, 5, 4, 4, 4, 4 };

  //交点の座標を最も原点に近い側のセル中心を原点とした座標で表す
  std::set< vector3<int> > tmp;
  for (int i = 0; i < 8; i++)
  {
    if ((cuts[i] & masks[i]) != 0LL)
    {
      if (ensCut(cuts[i], xdirs[i]) == 1)
      {
        int x = getBit9(cuts[i], xdirs[i]);
        if (x == 0)
        {
          for (int l = 0; l < 3; l++)
          {
            data[i][l] = 0;
          }
        }
        if (xdirs[i] % 2 == 0)
        {
          x = 511 - x;
        }
        tmp.insert(vector3<int>(x, (ydirs[i] + 1) % 2 * 511, (zdirs[i] + 1) % 2 * 511));
      }
      if (ensCut(cuts[i], ydirs[i]) == 1)
      {
        int y = getBit9(cuts[i], ydirs[i]);
        if (y == 0)
        {
          for (int l = 0; l < 3; l++)
          {
            data[i][l] = 0;
          }
        }
        if (ydirs[i] % 2 == 0)
        {
          y = 511 - y;
        }
        tmp.insert(vector3<int>((xdirs[i] + 1) % 2 * 511, y, (zdirs[i] + 1) % 2 * 511));
      }
      if (ensCut(cuts[i], zdirs[i]) == 1)
      {
        int z = getBit9(cuts[i], zdirs[i]);
        if (z == 0)
        {
          for (int l = 0; l < 3; l++)
          {
            data[i][l] = 0;
          }
        }
        if (zdirs[i] % 2 == 0)
        {
          z = 511 - z;
        }
        tmp.insert(vector3<int>((xdirs[i] + 1) % 2 * 511, (ydirs[i] + 1) % 2 * 511, z));
      }
    }
  }

  std::vector< vector3<int> > coords(tmp.begin(), tmp.end());

  //交点が1つだけの場合、交点は必ずどこかのセル中心と一致しなければならない。
  if (coords.size() == 1)
  {
    if (
        (coords[0].x != 0 && coords[0].x != 511) ||
        (coords[0].y != 0 && coords[0].y != 511) ||
        (coords[0].z != 0 && coords[0].z != 511))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("invalid cut plane 1");
    }
    //交点と重なっているセルの流速は0になっているので、通常通りトリリニア補間を呼ぶ
    return TrilinearInterpolate(dval, data, ip, jp, kp, im, jm, km);
  }

  //交点が2つだけの場合、交点を結ぶ直線はどこかのエッジと一致しなければならない
  if (coords.size() == 2)
  {
    if (
        (coords[0].x != coords[1].x) &&
        (coords[0].y != coords[1].y) &&
        (coords[0].z != coords[1].z))
    {
      LPT::LPT_LOG::GetInstance()->ERROR("invalid cut plane 2");
    }
    //交点と重なっているセルの流速は0になっているので、通常通りトリリニア補間を呼ぶ
    return TrilinearInterpolate(dval, data, ip, jp, kp, im, jm, km);
  }

  //交点が3点以上の場合は、最初の3点を使ってcut面の法線ベクトルを作成する
  long long normal_vector[3];
  if (coords.size() >= 3)
  {
    long long v0[3] = { coords[0].x, coords[0].x, coords[0].x };
    long long v1[3] = { coords[1].y, coords[1].y, coords[1].y };
    long long v2[3] = { coords[2].z, coords[2].z, coords[2].z };

    Utility::xmy< 3 >(v0, v1);
    Utility::xmy< 3 >(v0, v2);
    Utility::cross_product(v1, v2, normal_vector);
    //交点が4点以上の時は4点目以降が同一平面にあるかチェック
    if (coords.size() >= 4)
    {
      bool debug_print_done = false;
      for (int i = 4; i < coords.size(); i++)
      {
        long long v[3] = { coords[i].x, coords[i].y, coords[i].z };
        Utility::xmy< 3 >(v0, v);
        long long tmp = Utility::dot< 3 >(normal_vector, v);
        if (Utility::dot< 3 >(normal_vector, v) != 0)
        {
          LPT::LPT_LOG::GetInstance()->ERROR("invalid cut plane 3");
          if (!debug_print_done)
          {
            debug_print_done = true;
            LPT::LPT_LOG::GetInstance()->ERROR("origin = ",Origin,3);
            LPT::LPT_LOG::GetInstance()->ERROR("pitch  = ",Pitch,3);
            for (std::vector< vector3<int> >::iterator it = coords.begin(); it != coords.end(); ++it)
            {
              std::stringstream ss;
              ss << "cut coordinate(quantitize) = " << it->x << "," << it->y << "," << it->z;
              LPT::LPT_LOG::GetInstance()->ERROR(ss.str());
              std::stringstream ss2;
              ss2 << "cut coord = " << Origin[0] + Pitch[0]*(im + double((it->x))/512) << ","
                                   << Origin[1] + Pitch[1]*(jm + double((it->y))/512) << "," 
                                   << Origin[2] + Pitch[2]*(km + double((it->z))/512)  ;
              LPT::LPT_LOG::GetInstance()->ERROR(ss2.str());
            }
            for (int i = 0; i < 8; i++)
            {
              std::stringstream ss;
              ss << static_cast< std::bitset< 64 > >(cuts[i]);
              LPT::LPT_LOG::GetInstance()->ERROR(ss.str());
            }
          }
          //本来はここで何か対策しないといけないけど、とりあえず計算を進めるために
          //通常の線形補間をしてみる。
          return TrilinearInterpolate(dval, data, ip, jp, kp, im, jm, km);
        }
      }
    }

    // 流速の補間に使うセルを決める
    // pは1つ目の交点座標から粒子座標へのベクトル
    long long p[3];
    p[0] = quantize9(ip) - v0[0];
    p[1] = quantize9(jp) - v0[1];
    p[2] = quantize9(kp) - v0[2];

    long long pn = Utility::dot< 3 >(p, normal_vector);
    // 粒子が壁面にめり込んでいる時は流速は0とする
    if (pn == 0)
    {
      for (int l = 0; l < 3; l++)
      {
        dval[l] = 0;
      }
      return true;
    }

    // 補間に用いるセルを決める
    for (int i = 0; i < 8; i++)
    {
      // cellは1つ目の交点座標から各セル中心へのベクトル
      long long cell[3];
      cell[0]       = (xdirs[i] + 1) % 2 * 511 - v0[0];
      cell[1]       = (ydirs[i] + 1) % 2 * 511 - v0[1];
      cell[2]       = (zdirs[i] + 1) % 2 * 511 - v0[2];
      long long tmp = Utility::dot< 3 >(cell, normal_vector);

      //cut面を挟んで反対側にあるセル中心の流速は0にする
      //cut面上の1点(1つ目の交点)から各セルおよび粒子へのベクトルとcut面の法線ベクトルの内積が
      //異符号だったら、cut面を挟んで反対側にあるはず
      if (tmp * pn <= 0)
      {
        for (int l = 0; l < 3; l++)
        {
          data[i][l] = 0;
        }
      }
    }

    // 各セル中心の係数(alpha)を決める
    int p2[3];
    p2[0] = quantize9(ip);
    p2[1] = quantize9(jp);
    p2[2] = quantize9(kp);
    REAL_TYPE alpha[8];

    /*
    plus
    alpha[cell] = ensCut(cuts[cell], xdirs[cell]) == 1? REAL_TYPE(getBit9(cuts[cell], xdirs[cell])-p2[0])/REAL_TYPE(getBit9(cuts[cell], xdirs[cell])):im
    alpha[cell]*= ensCut(cuts[cell], ydirs[cell]) == 1? REAL_TYPE(getBit9(cuts[cell], ydirs[cell])-p2[1])/REAL_TYPE(getBit9(cuts[cell], ydirs[cell])):jm
    alpha[cell]*= ensCut(cuts[cell], zdirs[cell]) == 1? REAL_TYPE(getBit9(cuts[cell], zdirs[cell])-p2[2])/REAL_TYPE(getBit9(cuts[cell], zdirs[cell])):km
    minus
    alpha[cell] = ensCut(cuts[cell], xdirs[cell]) == 1? REAL_TYPE(p2[0]-(512-getBit9(cuts[cell], xdirs[cell])))/REAL_TYPE(getBit9(cuts[cell], xdirs[cell])):ip
    alpha[cell]*= ensCut(cuts[cell], ydirs[cell]) == 1? REAL_TYPE(p2[1]-(512-getBit9(cuts[cell], ydirs[cell])))/REAL_TYPE(getBit9(cuts[cell], ydirs[cell])):jp
    alpha[cell]*= ensCut(cuts[cell], zdirs[cell]) == 1? REAL_TYPE(p2[2]-(512-getBit9(cuts[cell], zdirs[cell])))/REAL_TYPE(getBit9(cuts[cell], zdirs[cell])):kp
    */

    //cell 0
    alpha[0] = ensCut(cuts[0], xdirs[0]) == 1 ? REAL_TYPE(getBit9(cuts[0], xdirs[0]) - p2[0]) / REAL_TYPE(getBit9(cuts[0], xdirs[0])) : im;
    alpha[0] *= ensCut(cuts[0], ydirs[0]) == 1 ? REAL_TYPE(getBit9(cuts[0], ydirs[0]) - p2[1]) / REAL_TYPE(getBit9(cuts[0], ydirs[0])) : jm;
    alpha[0] *= ensCut(cuts[0], zdirs[0]) == 1 ? REAL_TYPE(getBit9(cuts[0], zdirs[0]) - p2[2]) / REAL_TYPE(getBit9(cuts[0], zdirs[0])) : km;

    //cell 1
    alpha[1] = ensCut(cuts[1], xdirs[1]) == 1 ? REAL_TYPE(p2[0] - (512 - getBit9(cuts[1], xdirs[1]))) / REAL_TYPE(getBit9(cuts[1], xdirs[1])) : ip;
    alpha[1] *= ensCut(cuts[1], ydirs[1]) == 1 ? REAL_TYPE(getBit9(cuts[1], ydirs[1]) - p2[1]) / REAL_TYPE(getBit9(cuts[1], ydirs[1])) : jm;
    alpha[1] *= ensCut(cuts[1], zdirs[1]) == 1 ? REAL_TYPE(getBit9(cuts[1], zdirs[1]) - p2[2]) / REAL_TYPE(getBit9(cuts[1], zdirs[1])) : km;

    //cel 2
    alpha[2] = ensCut(cuts[2], xdirs[2]) == 1 ? REAL_TYPE(p2[0] - (512 - getBit9(cuts[2], xdirs[2]))) / REAL_TYPE(getBit9(cuts[2], xdirs[2])) : ip;
    alpha[2] *= ensCut(cuts[2], ydirs[2]) == 1 ? REAL_TYPE(p2[1] - (512 - getBit9(cuts[2], ydirs[2]))) / REAL_TYPE(getBit9(cuts[2], ydirs[2])) : jp;
    alpha[2] *= ensCut(cuts[2], zdirs[2]) == 1 ? REAL_TYPE(getBit9(cuts[2], zdirs[2]) - p2[2]) / REAL_TYPE(getBit9(cuts[2], zdirs[2])) : km;

    //cell 3
    alpha[3] = ensCut(cuts[3], xdirs[3]) == 1 ? REAL_TYPE(getBit9(cuts[3], xdirs[3]) - p2[0]) / REAL_TYPE(getBit9(cuts[3], xdirs[3])) : im;
    alpha[3] *= ensCut(cuts[3], ydirs[3]) == 1 ? REAL_TYPE(p2[1] - (512 - getBit9(cuts[3], ydirs[3]))) / REAL_TYPE(getBit9(cuts[3], ydirs[3])) : jp;
    alpha[3] *= ensCut(cuts[3], zdirs[3]) == 1 ? REAL_TYPE(getBit9(cuts[3], zdirs[3]) - p2[2]) / REAL_TYPE(getBit9(cuts[3], zdirs[3])) : km;

    //cell 4
    alpha[4] = ensCut(cuts[4], xdirs[4]) == 1 ? REAL_TYPE(getBit9(cuts[4], xdirs[4]) - p2[0]) / REAL_TYPE(getBit9(cuts[4], xdirs[4])) : im;
    alpha[4] *= ensCut(cuts[4], ydirs[4]) == 1 ? REAL_TYPE(getBit9(cuts[4], ydirs[4]) - p2[1]) / REAL_TYPE(getBit9(cuts[4], ydirs[4])) : jm;
    alpha[4] *= ensCut(cuts[4], zdirs[4]) == 1 ? REAL_TYPE(p2[2] - (512 - getBit9(cuts[4], zdirs[4]))) / REAL_TYPE(getBit9(cuts[4], zdirs[4])) : kp;

    //cell 5
    alpha[5] = ensCut(cuts[5], xdirs[5]) == 1 ? REAL_TYPE(p2[0] - (512 - getBit9(cuts[5], xdirs[5]))) / REAL_TYPE(getBit9(cuts[5], xdirs[5])) : ip;
    alpha[5] *= ensCut(cuts[5], ydirs[5]) == 1 ? REAL_TYPE(getBit9(cuts[5], ydirs[5]) - p2[1]) / REAL_TYPE(getBit9(cuts[5], ydirs[5])) : jm;
    alpha[5] *= ensCut(cuts[5], zdirs[5]) == 1 ? REAL_TYPE(p2[2] - (512 - getBit9(cuts[5], zdirs[5]))) / REAL_TYPE(getBit9(cuts[5], zdirs[5])) : kp;

    //cell 6
    alpha[6] =  ensCut(cuts[6], xdirs[6]) == 1 ? REAL_TYPE(p2[0] - (512 - getBit9(cuts[6], xdirs[6]))) / REAL_TYPE(getBit9(cuts[6], xdirs[6])) : ip;
    alpha[6] *= ensCut(cuts[6], ydirs[6]) == 1 ? REAL_TYPE(p2[1] - (512 - getBit9(cuts[6], ydirs[6]))) / REAL_TYPE(getBit9(cuts[6], ydirs[6])) : jp;
    alpha[6] *= ensCut(cuts[6], zdirs[6]) == 1 ? REAL_TYPE(p2[2] - (512 - getBit9(cuts[6], zdirs[6]))) / REAL_TYPE(getBit9(cuts[6], zdirs[6])) : kp;

    //cell 7
    alpha[7] = ensCut(cuts[7], xdirs[7]) == 1 ? REAL_TYPE(getBit9(cuts[7], xdirs[7]) - p2[0]) / REAL_TYPE(getBit9(cuts[7], xdirs[7])) : im;
    alpha[7] *= ensCut(cuts[7], ydirs[7]) == 1 ? REAL_TYPE(p2[1] - (512 - getBit9(cuts[7], ydirs[7]))) / REAL_TYPE(getBit9(cuts[7], ydirs[7])) : jp;
    alpha[7] *= ensCut(cuts[7], zdirs[7]) == 1 ? REAL_TYPE(p2[2] - (512 - getBit9(cuts[7], zdirs[7]))) / REAL_TYPE(getBit9(cuts[7], zdirs[7])) : kp;

    for (int l = 0; l < 3; l++)
    {
      dval[l] = 0;
    }
    //係数をかけながら周辺セルの流速を足し込む
    for (int i = 0; i < 8; i++)
    {
      for (int l = 0; l < 3; l++)
      {
        dval[l] += alpha[i] * data[i][l];
      }
    }
  }

  return true;
}

#define INDEX(i, j, k, l) ((i) + (j)*DataBlock.BlockSize[0] + (k)*DataBlock.BlockSize[0] * DataBlock.BlockSize[1] + (l)*DataBlock.BlockSize[0] * DataBlock.BlockSize[1] * DataBlock.BlockSize[2])
bool Interpolator::InterpolateData(const DSlib::DataBlock& DataBlock, const REAL_TYPE x_I[3], REAL_TYPE dval[3], const std::string id)
{
  if (!DataBlock.Data) return false;
  int i = int(x_I[0]);
  int j = int(x_I[1]);
  int k = int(x_I[2]);

  REAL_TYPE ip = x_I[0] - (REAL_TYPE)i;
  REAL_TYPE jp = x_I[1] - (REAL_TYPE)j;
  REAL_TYPE kp = x_I[2] - (REAL_TYPE)k;
  REAL_TYPE im = (REAL_TYPE)(i + 1) - x_I[0];
  REAL_TYPE jm = (REAL_TYPE)(j + 1) - x_I[1];
  REAL_TYPE km = (REAL_TYPE)(k + 1) - x_I[2];

  REAL_TYPE data[8][3];

  for (int l = 0; l < 3; l++)
  {
    data[0][l] = DataBlock.Data[INDEX(i, j, k, l)];
    data[1][l] = DataBlock.Data[INDEX(i + 1, j, k, l)];
    data[2][l] = DataBlock.Data[INDEX(i + 1, j + 1, k, l)];
    data[3][l] = DataBlock.Data[INDEX(i, j + 1, k, l)];
    data[4][l] = DataBlock.Data[INDEX(i, j, k + 1, l)];
    data[5][l] = DataBlock.Data[INDEX(i + 1, j, k + 1, l)];
    data[6][l] = DataBlock.Data[INDEX(i + 1, j + 1, k + 1, l)];
    data[7][l] = DataBlock.Data[INDEX(i, j + 1, k + 1, l)];
  }

  if (DataBlock.has_cut)
  {
    const long long cuts[8] = { DataBlock.d_cut[INDEX(i, j, k, 0)],
                                DataBlock.d_cut[INDEX(i + 1, j, k, 0)],
                                DataBlock.d_cut[INDEX(i + 1, j + 1, k, 0)],
                                DataBlock.d_cut[INDEX(i, j + 1, k, 0)],
                                DataBlock.d_cut[INDEX(i, j, k + 1, 0)],
                                DataBlock.d_cut[INDEX(i + 1, j, k + 1, 0)],
                                DataBlock.d_cut[INDEX(i + 1, j + 1, k + 1, 0)],
                                DataBlock.d_cut[INDEX(i, j + 1, k + 1, 0)] };

    long long has_cut_local = (cuts[0] & masks[0]) |
                              (cuts[1] & masks[1]) |
                              (cuts[2] & masks[2]) |
                              (cuts[3] & masks[3]) |
                              (cuts[4] & masks[4]) |
                              (cuts[5] & masks[5]) |
                              (cuts[6] & masks[6]) |
                              (cuts[7] & masks[7]);

    if (has_cut_local != 0)
    {
      LPT::LPT_LOG::GetInstance()->LOG("interpolate with cut ", id);
      return InterpolateWithCut(dval, data, cuts, ip, jp, kp, im, jm, km, DataBlock.Origin, DataBlock.Pitch);
    }
  }

  TrilinearInterpolate(dval, data, ip, jp, kp, im, jm, km);
  return true;
#undef INDEX
}

void Interpolator::ConvXtoI(const REAL_TYPE x_g[3], REAL_TYPE x_l[3], const REAL_TYPE orig[3], const REAL_TYPE pitch[3])
{
  static const int halo = DSlib::DecompositionManager::GetInstance()->GetGuideCellSize();
  x_l[0]                = ((x_g[0] - orig[0]) / pitch[0] + (halo - 0.5));
  x_l[1]                = ((x_g[1] - orig[1]) / pitch[1] + (halo - 0.5));
  x_l[2]                = ((x_g[2] - orig[2]) / pitch[2] + (halo - 0.5));
}

void Interpolator::ConvItoX(const REAL_TYPE x_l[3], REAL_TYPE x_g[3], const REAL_TYPE orig[3], const REAL_TYPE pitch[3])
{
  static const int halo = DSlib::DecompositionManager::GetInstance()->GetGuideCellSize();
  x_g[0]                = orig[0] + (x_l[0] - (halo - 0.5)) * pitch[0];
  x_g[1]                = orig[1] + (x_l[1] - (halo - 0.5)) * pitch[1];
  x_g[2]                = orig[2] + (x_l[2] - (halo - 0.5)) * pitch[2];
}
} // namespace PPlib
