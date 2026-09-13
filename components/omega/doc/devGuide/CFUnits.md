(omega-dev-cfunits)=

## CFUnits

Omega fields carry their units as a string attribute in the plain form used
by the CF standard name table: unit symbols separated by single spaces, each
with an optional signed integer exponent (`m s-1`, `kg m-3`, `W m-2`), and `1`
for a dimensionless quantity. When a field is derived from others, its units
must be derived too. The `CFUnits` class parses such a string, supports the
product, quotient and integer power of units expressions, and formats the
result back in the plain form.

The string interfaces cover the common cases:
```c++
#include "CFUnits.h"

std::string Flux   = CFUnits::multiply("m s-1", "m2");   // "m3 s-1"
std::string Speed  = CFUnits::divide("m3 s-1", "m2");    // "m s-1"
std::string Energy = CFUnits::power("m s-1", 2);         // "m2 s-2"
```
Exponents of the same symbol are merged (`m s-1` times `m` is `m2 s-1`) and
a symbol whose exponent reaches zero is dropped, so `kg m-3` times `m3` is
`kg`. In the result, symbols with positive exponents come first, then those
with negative exponents, each in the order they first appeared.

A field with no units attribute has *unknown* units, represented by the empty
string. Unknown units propagate through every operation, so a field derived
from one without units has no units either, rather than wrong ones.

The parser also accepts the udunits caret and slash spellings (`m^2`, `m/s`),
which some older field definitions still use, but always emits the plain
form. Any other string is an error: the string interfaces abort, and the
value interface reports it so a caller can decide:
```c++
CFUnits Units;
Error Err = CFUnits::parse("m s-1", Units);  // Fail for a malformed string
CFUnits Area;
CFUnits::parse("m2", Area);
std::string Flux = (Units * Area).str();     // "m3 s-1"
```
`CFUnits` also provides `operator/`, `pow(int)`, `isUnknown()`,
`isDimensionless()` and an order-insensitive `operator==`.

A malformed units string is a programming error in a field definition, and
aborting when the derived field is created stops the run at initialization
rather than writing output with a wrong or missing attribute.
