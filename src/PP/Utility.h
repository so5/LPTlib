/*
 * LPTlib
 * Lagrangian Particle Tracking library
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef PPLIB_UTILITY_H
#define PPLIB_UTILITY_H
#include <cmath>
#include <algorithm>
#include "LPT_LogOutput.h"

namespace PPlib
{
namespace Utility
{
//! @brief 与えられたNを素因数分解して、結果をFactorに収める
//! @param N      [in]  素因数分解を行なう数
//! @param Factor [out] 結果を格納するコンテナ
inline void Factorize(const int& argN, std::vector< std::pair< int, int > >* Factor)
{
  //念のためFactorを空にする
  Factor->clear();

  //Nが0以下なら何もしない
  if (argN <= 0)
  {
    return;
  }

  //Nが正の数ならまず1を代入
  Factor->push_back(std::pair< int, int >(1, 1));

  //Nが1の時は終了
  if (argN == 1)
  {
    return;
  }

  int N = argN;

  //2で何回割れるか？
  int pow2 = 0;
  while (N % 2 == 0)
  {
    N /= 2;
    pow2++;
  }
  if (pow2 > 0)
  {
    Factor->push_back(std::pair< int, int >(2, pow2));
  }
  if (N == 1)
  {
    return;
  }

  //3以上sqrt(N)以下の奇数の因数を探す
  for (int i = 3; i <= std::sqrt((double)N); i += 2)
  {
    int pow = 0;

    while (N % i == 0)
    {
      N /= i;
      pow++;
    }
    if (pow > 0)
    {
      Factor->push_back(std::pair< int, int >(i, pow));
    }
    if (N == 1)
    {
      return;
    }
  }

  //余りがあれば追加して終了
  Factor->push_back(std::pair< int, int >(N, 1));
}

//! @brief N^Mを計算して返す
//
//memo: cmath内のpow関数は浮動小数点数用のみしか無い
inline int pow(const int& N, const int& M)
{
  int rt = 1;
  for (int i = 0; i < M; i++)
    rt *= N;
  return rt;
}

//! @brief 与えられたNを2^p0*3^p1*5^p2*Reminder の形に分解し、2,3,5のべき数と余りを返す
//! @param N   [in]  素因数分解を行なう数
//! @param Pow [out] 結果を格納する配列
inline void Factorize235(const int& N, int Pow[4])
{
  //Pow[]を初期化
  Pow[0] = 0;
  Pow[1] = 0;
  Pow[2] = 0;
  Pow[3] = N;

  //2で何回割れるか？
  while (Pow[3] % 2 == 0)
  {
    Pow[3] /= 2;
    Pow[0] += 1;
  }

  //3で何回割れるか？
  while (Pow[3] % 3 == 0)
  {
    Pow[3] /= 3;
    Pow[1] += 1;
  }

  //5で何回割れるか？
  while (Pow[3] % 5 == 0)
  {
    Pow[3] /= 5;
    Pow[2] += 1;
  }

  //Pow[3]は余り
}

//! point1とpoint2をnum_points分割した値を昇順にソートしてrtに格納する
inline void DivideLine1D(std::vector< REAL_TYPE >* rt, const int& num_points, const REAL_TYPE& point1, const REAL_TYPE& point2)
{
  if (num_points == 1)
  {
    rt->push_back(point1);
  } else {
    for (int i = 0; i < num_points; i++)
    {
      if (point1 != point2)
      {
        rt->push_back(point1 + (point2 - point1) / (num_points - 1) * i);
      } else {
        rt->push_back(point1);
      }
    }
  }
  std::sort(rt->begin(), rt->end());
}

//! NxMxKの領域をNBxMBxKB に分割する
//
//余りは許容し、NBxMBxKB<MaxPoints となるように適当に調整する
inline void DetermineBlockSize(int* arg_NB, int* arg_MB, int* arg_KB, const int& MaxPoints, const int& N, const int& M, const int& K)
{
  int& NB = *arg_NB;
  int& MB = *arg_MB;
  int& KB = *arg_KB;

  std::vector< std::pair< int, int > > Factor;
  Factorize(MaxPoints, &Factor);
  NB = 1;
  MB = 1;
  KB = 1;
  for (std::vector< std::pair< int, int > >::iterator it = Factor.begin(); it != Factor.end(); ++it)
  {
    NB *= pow(it->first, (it->second) / 3);
    MB *= pow(it->first, (it->second) / 3);
    KB *= pow(it->first, (it->second) / 3);
    if ((it->second) % 3 != 0)
    {
      if (NB <= MB && NB <= KB)
      {
        NB *= pow(it->first, (it->second) % 3);
      } else if (MB <= KB) {
        MB *= pow(it->first, (it->second) % 3);
      } else {
        KB *= pow(it->first, (it->second) % 3);
      }
    }
  }

  LPT::LPT_LOG::GetInstance()->LOG("initial NB = ", NB);
  LPT::LPT_LOG::GetInstance()->LOG("initial MB = ", MB);
  LPT::LPT_LOG::GetInstance()->LOG("initial KB = ", KB);

  // NB,MB,KBがそれぞれN,M,Kを越えてたらN,M,Kに置き換える
  if (NB > N) NB = N;
  if (MB > M) MB = M;
  if (KB > K) KB = K;

  //もしNB, MB, KBを大きくしても条件を満たせるなら大きくする
  while ((NB + 1) * MB * KB <= MaxPoints && NB + 1 <= N)
    ++NB;
  while (NB * (MB + 1) * KB <= MaxPoints && MB + 1 <= M)
    ++MB;
  while (NB * MB * (KB + 1) <= MaxPoints && KB + 1 <= K)
    ++KB;

  LPT::LPT_LOG::GetInstance()->LOG("NB = ", NB);
  LPT::LPT_LOG::GetInstance()->LOG("MB = ", MB);
  LPT::LPT_LOG::GetInstance()->LOG("KB = ", KB);
}

//! 引数で渡された3次元ベクトルを単位ベクトルに変換する
template < typename T >
inline void NormalizeVector(T* v)
{
  // オーバーフロー対策としていずれかの要素が2^15を越えていた場合は1回その値で割っておく
  for (int i = 0; i < 3; i++)
  {
    if (v[i] > 32768)
    {
      v[0] /= v[i];
      v[1] /= v[i];
      v[2] /= v[i];
    }
  }
  const T length = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (length != 0)
  {
    v[0] = v[0] / length;
    v[1] = v[1] / length;
    v[2] = v[2] / length;
  }
}

//! 引数で渡されたベクトルの内積を計算する
template < size_t N, typename T >
inline T dot(T* v1, T* v2)
{
  T rt;
  for (size_t i = 0; i < N; i++)
  {
    rt += v1[i] * v2[i];
  }
  return rt;
}

//! 引数で渡された3次元ベクトルの外積を計算する
template < typename T >
inline void cross_product(T* v1, T* v2, T* cross_product)
{
  cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

//! 二つのベクトルの差分を取り2つ目の配列に格納する
//
//BLASの*axpyをA=-1, incx=1, incy=1として呼び出したのと同等の処理
template < size_t N, typename T >
inline void xmy(T* X, T* Y)
{
  for (size_t i = 0; i < N; i++)
  {
    Y[i] -= X[i];
  }
}
}
} // namespace PPlib
#endif
