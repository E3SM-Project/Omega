#ifndef OMEGA_EOSCONSTANTS_H
#define OMEGA_EOSCONSTANTS_H
//===-- ocn/Eos.h - Equation of State --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating specific volume
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the specific volume based on the choice of EOS
//
//===----------------------------------------------------------------------===//

namespace OMEGA {

#define GSW_SPECVOL_COEFFICIENTS
const double V000 = 1.0769995862e-03;
const double V100 = -3.1038981976e-04;
const double V200 = 6.6928067038e-04;
const double V300 = -8.5047933937e-04;
const double V400 = 5.8086069943e-04;
const double V500 = -2.1092370507e-04;
const double V600 = 3.1932457305e-05;
const double V010 = -1.5649734675e-05;
const double V110 = 3.5009599764e-05;
const double V210 = -4.3592678561e-05;
const double V310 = 3.4532461828e-05;
const double V410 = -1.1959409788e-05;
const double V510 = 1.3864594581e-06;
const double V020 = 2.7762106484e-05;
const double V120 = -3.7435842344e-05;
const double V220 = 3.5907822760e-05;
const double V320 = -1.8698584187e-05;
const double V420 = 3.8595339244e-06;
const double V030 = -1.6521159259e-05;
const double V130 = 2.4141479483e-05;
const double V230 = -1.4353633048e-05;
const double V330 = 2.2863324556e-06;
const double V040 = 6.9111322702e-06;
const double V140 = -8.7595873154e-06;
const double V240 = 4.3703680598e-06;
const double V050 = -8.0539615540e-07;
const double V150 = -3.3052758900e-07;
const double V060 = 2.0543094268e-07;
const double V001 = -1.6784136540e-05;
const double V101 = 2.4262468747e-05;
const double V201 = -3.4792460974e-05;
const double V301 = 3.7470777305e-05;
const double V401 = -1.7322218612e-05;
const double V501 = 3.0927427253e-06;
const double V011 = 1.8505765429e-05;
const double V111 = -9.5677088156e-06;
const double V211 = 1.1100834765e-05;
const double V311 = -9.8447117844e-06;
const double V411 = 2.5909225260e-06;
const double V021 = -1.1716606853e-05;
const double V121 = -2.3678308361e-07;
const double V221 = 2.9283346295e-06;
const double V321 = -4.8826139200e-07;
const double V031 = 7.9279656173e-06;
const double V131 = -3.4558773655e-06;
const double V231 = 3.1655306078e-07;
const double V041 = -3.4102187482e-06;
const double V141 = 1.2956717783e-06;
const double V051 = 5.0736766814e-07;
const double V002 = 3.0623833435e-06;
const double V102 = -5.8484432984e-07;
const double V202 = -4.8122251597e-06;
const double V302 = 4.9263106998e-06;
const double V402 = -1.7811974727e-06;
const double V012 = -1.1736386731e-06;
const double V112 = -5.5699154557e-06;
const double V212 = 5.4620748834e-06;
const double V312 = -1.3544185627e-06;
const double V022 = 2.1305028740e-06;
const double V122 = 3.9137387080e-07;
const double V222 = -6.5731104067e-07;
const double V032 = -4.6132540037e-07;
const double V132 = 7.7618888092e-09;
const double V042 = -6.3352916514e-08;
const double V003 = -3.8088938393e-07;
const double V103 = 3.6310188515e-07;
const double V203 = 1.6746303780e-08;
const double V013 = -3.6527006553e-07;
const double V113 = -2.7295696237e-07;
const double V023 = 2.8695905159e-07;
const double V004 = 8.8302421514e-08;
const double V104 = -1.1147125423e-07;
const double V014 = 3.1454099902e-07;
const double V005 = 4.2369007180e-09;

} // namespace OMEGA
#endif
