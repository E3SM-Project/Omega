//===-- Test driver for OMEGA CFUnits ---------------------------*- C++ -*-===//
//
/// \file
/// \brief Test driver for the CFUnits units algebra
///
/// Checks that units strings in the CF plain form and the udunits caret and
/// slash spellings parse and format as expected, that products, quotients
/// and powers give the plain form with merged exponents, that unknown
/// (empty) units propagate, and that malformed strings are rejected.
//
//===----------------------------------------------------------------------===//

#include "CFUnits.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"

#include <string>

using namespace OMEGA;

namespace {

int NumFailed = 0;

// Checks that a units string parses and formats to the expected plain form
void checkParse(const std::string &Input, const std::string &Expected) {
   CFUnits Units;
   Error Err = CFUnits::parse(Input, Units);
   if (!Err.isSuccess()) {
      LOG_ERROR("CFUnitsTest: '{}' failed to parse: {}", Input, Err.Msg);
      ++NumFailed;
   } else if (Units.str() != Expected) {
      LOG_ERROR("CFUnitsTest: '{}' formatted as '{}', expected '{}'", Input,
                Units.str(), Expected);
      ++NumFailed;
   }
}

// Checks that a malformed units string is rejected
void checkInvalid(const std::string &Input) {
   CFUnits Units;
   Error Err = CFUnits::parse(Input, Units);
   if (!Err.isFail()) {
      LOG_ERROR("CFUnitsTest: '{}' parsed as '{}' but should be rejected",
                Input, Units.str());
      ++NumFailed;
   } else if (!Units.isUnknown()) {
      LOG_ERROR("CFUnitsTest: '{}' was rejected but left units '{}'", Input,
                Units.str());
      ++NumFailed;
   }
}

// Checks a string result against its expected value
void checkResult(const std::string &Label, const std::string &Result,
                 const std::string &Expected) {
   if (Result != Expected) {
      LOG_ERROR("CFUnitsTest: {} gave '{}', expected '{}'", Label, Result,
                Expected);
      ++NumFailed;
   }
}

// Checks a boolean condition
void checkTrue(const std::string &Label, bool Condition) {
   if (!Condition) {
      LOG_ERROR("CFUnitsTest: {} is false", Label);
      ++NumFailed;
   }
}

} // namespace

int main(int argc, char **argv) {

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   initLogging(DefEnv);

   LOG_INFO("----- CFUnitsTest -----");

   // Parsing and formatting in the CF plain form
   checkParse("", "");
   checkParse("   ", "");
   checkParse("1", "1");
   checkParse("m", "m");
   checkParse(" m ", "m");
   checkParse("m s-1", "m s-1");
   checkParse("kg m-3", "kg m-3");
   checkParse("m3 kg-1 K-1", "m3 kg-1 K-1");
   checkParse("degree_C", "degree_C");
   checkParse("m s-1 m", "m2 s-1");
   checkParse("m m-1", "1");
   checkParse("m1", "m");
   checkParse("m0 s", "s");

   // Positive exponents are written before negative ones
   checkParse("s-1 m", "m s-1");
   checkParse("m-3 kg", "kg m-3");

   // The udunits caret and slash spellings are accepted
   checkParse("m^2", "m2");
   checkParse("m^-1 s^-1", "m-1 s-1");
   checkParse("m/s", "m s-1");
   checkParse("m/s/s", "m s-2");
   checkParse("1/s", "s-1");
   checkParse("kg m-3/s", "kg m-3 s-1");
   checkParse("m^2/s", "m2 s-1");

   // Malformed strings are rejected and leave unknown units
   checkInvalid("N m^{-2}");
   checkInvalid("m*s");
   checkInvalid("m.s");
   checkInvalid("m^");
   checkInvalid("m^+2");
   checkInvalid("m--1");
   checkInvalid("m2s");
   checkInvalid("2 m");
   checkInvalid("1e6 m3 s-1");
   checkInvalid("m/");
   checkInvalid("/s");
   checkInvalid("m//s");
   checkInvalid("(m s)");
   checkInvalid("m,s");

   // Queries
   {
      CFUnits Unknown, Dimless, Meters;
      CFUnits::parse("", Unknown);
      CFUnits::parse("1", Dimless);
      CFUnits::parse("m", Meters);
      checkTrue("empty is unknown", Unknown.isUnknown());
      checkTrue("empty is not dimensionless", !Unknown.isDimensionless());
      checkTrue("1 is dimensionless", Dimless.isDimensionless());
      checkTrue("1 is not unknown", !Dimless.isUnknown());
      checkTrue("m is not unknown", !Meters.isUnknown());
      checkTrue("m is not dimensionless", !Meters.isDimensionless());
      checkTrue("default is unknown", CFUnits().isUnknown());
   }

   // Equality ignores the order of factors and distinguishes unknown
   {
      CFUnits A, B, C, Unknown;
      CFUnits::parse("m s-1", A);
      CFUnits::parse("s-1 m", B);
      CFUnits::parse("m s-2", C);
      checkTrue("m s-1 == s-1 m", A == B);
      checkTrue("m s-1 != m s-2", A != C);
      checkTrue("m s-1 != unknown", A != Unknown);
      checkTrue("unknown == unknown", Unknown == CFUnits());
   }

   // Products, quotients and powers of strings
   checkResult("m s-1 * m2", CFUnits::multiply("m s-1", "m2"), "m3 s-1");
   checkResult("kg m-3 * m3", CFUnits::multiply("kg m-3", "m3"), "kg");
   checkResult("m s-1 * kg m-3", CFUnits::multiply("m s-1", "kg m-3"),
               "kg m-2 s-1");
   checkResult("m * m-1", CFUnits::multiply("m", "m-1"), "1");
   checkResult("1 * m", CFUnits::multiply("1", "m"), "m");
   checkResult("degree_C * m2", CFUnits::multiply("degree_C", "m2"),
               "degree_C m2");
   checkResult("W m-2 * m2", CFUnits::multiply("W m-2", "m2"), "W");
   checkResult("m3 s-1 / m2", CFUnits::divide("m3 s-1", "m2"), "m s-1");
   checkResult("m / m", CFUnits::divide("m", "m"), "1");
   checkResult("1 / s", CFUnits::divide("1", "s"), "s-1");
   checkResult("(m s-1)^2", CFUnits::power("m s-1", 2), "m2 s-2");
   checkResult("(m s-1)^-1", CFUnits::power("m s-1", -1), "s m-1");
   checkResult("m^0", CFUnits::power("m", 0), "1");

   // Unknown units propagate through every operation
   checkResult("unknown * m2", CFUnits::multiply("", "m2"), "");
   checkResult("m2 * unknown", CFUnits::multiply("m2", ""), "");
   checkResult("unknown / m2", CFUnits::divide("", "m2"), "");
   checkResult("m2 / unknown", CFUnits::divide("m2", ""), "");
   checkResult("unknown^3", CFUnits::power("", 3), "");
   checkResult("unknown^0", CFUnits::power("", 0), "");

   // The same algebra on parsed values
   {
      CFUnits Vel, Area, Dens;
      CFUnits::parse("m s-1", Vel);
      CFUnits::parse("m2", Area);
      CFUnits::parse("kg m-3", Dens);
      checkResult("value product", (Vel * Area).str(), "m3 s-1");
      checkResult("value quotient", (Vel / Area).str(), "m-1 s-1");
      checkResult("value chain", (Vel * Area * Dens).str(), "kg s-1");
      checkResult("value power", Dens.pow(-1).str(), "m3 kg-1");
      checkResult("value unknown", (Vel * CFUnits()).str(), "");
   }

   if (NumFailed > 0) {
      ABORT_ERROR("CFUnitsTest: {} checks failed", NumFailed);
   }

   LOG_INFO("----- CFUnitsTest Successful -----");

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return 0;
}
