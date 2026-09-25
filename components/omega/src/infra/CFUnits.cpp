//===-- infra/CFUnits.cpp - units algebra for CF units strings --*- C++ -*-===//
//
// Implementation of the CFUnits class. A units string is split on slashes
// into a numerator and denominators, each of those on whitespace into
// factors, and each factor into a unit symbol and an optional integer
// exponent. Factors are kept in order of first appearance so that the
// formatted result reads like the conventional spelling (kg m-3, not
// m-3 kg).
//
//===----------------------------------------------------------------------===//

#include "CFUnits.h"
#include "Logging.h"

#include <cctype>
#include <map>
#include <sstream>

namespace OMEGA {

//------------------------------------------------------------------------------
// Default constructor creates unknown units
CFUnits::CFUnits() : Known(false) {}

//------------------------------------------------------------------------------
// Parses a units string into its factors
Error CFUnits::parse(const std::string &UnitsStr, // [in] units string
                     CFUnits &Result // [out] parsed units expression
) {

   Error Err;

   Result = CFUnits();

   // Trim surrounding whitespace; an empty string is the unknown units of a
   // field with no units attribute
   auto First = UnitsStr.find_first_not_of(" \t");
   if (First == std::string::npos)
      return Err;
   auto Last           = UnitsStr.find_last_not_of(" \t");
   std::string Trimmed = UnitsStr.substr(First, Last - First + 1);

   Result.Known = true;

   // A trailing slash has no denominator after it; getline would skip it
   if (Trimmed.back() == '/') {
      Result = CFUnits();
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "CFUnits: cannot parse units string '{}': "
                   "empty denominator",
                   UnitsStr);
   }

   // Each slash starts a denominator: every factor after it divides
   int Sign = 1;
   std::string Segment;
   std::istringstream Segments(Trimmed);
   bool FirstSegment = true;
   while (std::getline(Segments, Segment, '/')) {

      // A slash with nothing before or after it is malformed
      std::istringstream Tokens(Segment);
      std::string Token;
      bool AnyToken = false;
      while (Tokens >> Token) {
         AnyToken = true;
         if (!Result.parseFactor(Token, Sign)) {
            Result = CFUnits();
            RETURN_ERROR(Err, ErrorCode::Fail,
                         "CFUnits: cannot parse units string '{}': "
                         "unexpected factor '{}'",
                         UnitsStr, Token);
         }
      }
      if (!AnyToken) {
         Result = CFUnits();
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "CFUnits: cannot parse units string '{}': "
                      "empty {}",
                      UnitsStr, FirstSegment ? "numerator" : "denominator");
      }

      FirstSegment = false;
      Sign         = -1;
   }

   return Err;

} // end parse

//------------------------------------------------------------------------------
// Parses one factor: a unit symbol of letters and underscores followed by an
// optional integer exponent, with or without a caret (m2, s-1, m^2, s^-1), or
// the digit 1 which contributes nothing
bool CFUnits::parseFactor(const std::string &Token, // [in] one factor
                          int Sign                  // [in] +1 or -1
) {

   if (Token == "1")
      return true;

   // Symbol: one or more letters or underscores
   std::size_t Pos = 0;
   while (Pos < Token.size() &&
          (std::isalpha(static_cast<unsigned char>(Token[Pos])) ||
           Token[Pos] == '_'))
      ++Pos;
   if (Pos == 0)
      return false;
   std::string Symbol = Token.substr(0, Pos);

   // Exponent: optional caret, optional minus sign, then digits to the end
   int Exponent = 1;
   if (Pos < Token.size()) {
      if (Token[Pos] == '^')
         ++Pos;
      int ExpSign = 1;
      if (Pos < Token.size() && Token[Pos] == '-') {
         ExpSign = -1;
         ++Pos;
      }
      if (Pos >= Token.size())
         return false;
      int Magnitude = 0;
      for (; Pos < Token.size(); ++Pos) {
         if (!std::isdigit(static_cast<unsigned char>(Token[Pos])))
            return false;
         Magnitude = 10 * Magnitude + (Token[Pos] - '0');
      }
      Exponent = ExpSign * Magnitude;
   }

   multiplyFactor(Symbol, Sign * Exponent);
   return true;

} // end parseFactor

//------------------------------------------------------------------------------
// Multiplies in one symbol raised to a power
void CFUnits::multiplyFactor(const std::string &Symbol, // [in] unit symbol
                             int Exponent               // [in] its power
) {

   if (Exponent == 0)
      return;

   for (auto Iter = Factors.begin(); Iter != Factors.end(); ++Iter) {
      if (Iter->first == Symbol) {
         Iter->second += Exponent;
         if (Iter->second == 0)
            Factors.erase(Iter);
         return;
      }
   }
   Factors.emplace_back(Symbol, Exponent);

} // end multiplyFactor

//------------------------------------------------------------------------------
// Parses a units string, aborting on failure
CFUnits CFUnits::parseOrAbort(const std::string &UnitsStr // [in] units string
) {
   CFUnits Result;
   Error Err = parse(UnitsStr, Result);
   CHECK_ERROR_ABORT(Err, "CFUnits: invalid units string");
   return Result;
}

//------------------------------------------------------------------------------
// String interfaces for the product, quotient and power
std::string CFUnits::multiply(const std::string &UnitsA, // [in] first units
                              const std::string &UnitsB  // [in] second units
) {
   return (parseOrAbort(UnitsA) * parseOrAbort(UnitsB)).str();
}

std::string CFUnits::divide(const std::string &UnitsA, // [in] numerator
                            const std::string &UnitsB  // [in] denominator
) {
   return (parseOrAbort(UnitsA) / parseOrAbort(UnitsB)).str();
}

std::string CFUnits::power(const std::string &UnitsA, // [in] units
                           int Exponent               // [in] power
) {
   return parseOrAbort(UnitsA).pow(Exponent).str();
}

//------------------------------------------------------------------------------
// Queries
bool CFUnits::isUnknown() const { return !Known; }

bool CFUnits::isDimensionless() const { return Known && Factors.empty(); }

//------------------------------------------------------------------------------
// Formats the units in the CF plain form
std::string CFUnits::str() const {

   if (!Known)
      return "";
   if (Factors.empty())
      return "1";

   std::string Result;
   for (int Pass = 0; Pass < 2; ++Pass) {
      for (const auto &Factor : Factors) {
         bool Positive = Factor.second > 0;
         if (Positive != (Pass == 0))
            continue;
         if (!Result.empty())
            Result += " ";
         Result += Factor.first;
         if (Factor.second != 1)
            Result += std::to_string(Factor.second);
      }
   }
   return Result;

} // end str

//------------------------------------------------------------------------------
// Algebra
CFUnits CFUnits::operator*(const CFUnits &Other) const {

   if (!Known || !Other.Known)
      return CFUnits();

   CFUnits Result = *this;
   for (const auto &Factor : Other.Factors)
      Result.multiplyFactor(Factor.first, Factor.second);
   return Result;
}

CFUnits CFUnits::operator/(const CFUnits &Other) const {
   return *this * Other.pow(-1);
}

CFUnits CFUnits::pow(int Exponent) const {

   if (!Known)
      return CFUnits();

   CFUnits Result;
   Result.Known = true;
   for (const auto &Factor : Factors)
      Result.multiplyFactor(Factor.first, Factor.second * Exponent);
   return Result;
}

//------------------------------------------------------------------------------
// Comparison, ignoring the order of factors
bool CFUnits::operator==(const CFUnits &Other) const {

   if (Known != Other.Known)
      return false;
   std::map<std::string, int> Mine(Factors.begin(), Factors.end());
   std::map<std::string, int> Theirs(Other.Factors.begin(),
                                     Other.Factors.end());
   return Mine == Theirs;
}

bool CFUnits::operator!=(const CFUnits &Other) const {
   return !(*this == Other);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
