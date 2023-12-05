//===-- Test driver for OMEGA Metadata -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Metadata
///
/// This driver tests the Metadata capabilities for the OMEGA
/// model
//
//===-----------------------------------------------------------------------===/

#include <iostream>
#include <limits>

#include "DataTypes.h"
#include "Logging.h"
#include "MetaData.h"

using namespace OMEGA;

void printResult(bool result, bool expected, std::string msg_pass,
                 std::string msg_fail) {

   if (result == expected) {
      std::cout << msg_pass << ": PASS" << std::endl;
   } else {
      std::cout << msg_fail << ": FAIL" << std::endl;
   }
}

void testMetaDim() {

   const std::string DimName{"MyDim"};
   const I4 DimValue{1};

   // test has()
   printResult(MetaDim::has(DimName), false, "'" + DimName + "' is not created",
               "'" + DimName + "' should not exist");

   // test create()
   auto Dim1 = MetaDim::create(DimName, DimValue);

   printResult(MetaDim::has(DimName), true, "'" + DimName + "' is created",
               "'" + DimName + "' should exist");

   // test get()
   auto Dim2 = MetaDim::get(DimName);

   printResult(Dim1 == Dim2, true, "get() returns correct instance.",
               "get() returns incorrect instance.");

   // test getLength()
   I4 Length;
   Dim1->getLength(Length);

   printResult(DimValue == Length, true, "getLength() returns correct length.",
               "getLength() returns incorrect length.");

   // test destroy()
   MetaDim::destroy(DimName);

   printResult(MetaDim::has(DimName), false,
               "'" + DimName + "' is destroyed correctly",
               "'" + DimName + "' is not destroyed");
}

void testMetaData() {

   const std::string ArrayName{"MyArray"};
   const std::string DimName{"MyDim"};
   const I4 DimValue{1};
   int ret;

   // test has()
   printResult(MetaData::has(CodeMeta), false,
               "'" + CodeMeta + "' is not created",
               "'" + CodeMeta + "' should not exist");

   // test create() 1
   auto Data1 = MetaData::create(CodeMeta);

   printResult(MetaData::has(CodeMeta), true, "'" + CodeMeta + "' is created",
               "'" + CodeMeta + "' should exist");

   // test create() 2
   auto Data2 = MetaData::create(SimMeta, {std::make_pair("Meta1", 1),
                                           std::make_pair("Meta2", 2),
                                           std::make_pair("Meta3", 3)});

   printResult(MetaData::has(SimMeta), true, "'" + SimMeta + "' is created",
               "'" + SimMeta + "' should exist");

   std::map<std::string, std::any> *VarMeta = Data2->getAllMetaData();
   std::string MetaName;
   int MetaVal, Count = 1;

   // use std::map functionality to loop through all metadata
   for (const auto &Pair : *VarMeta) {
      MetaName = Pair.first; // retrieve name part of meta pair
      MetaVal  = std::any_cast<int>(Pair.second); // retrieve value part

      printResult(MetaName == "Meta" + std::to_string(MetaVal), true,
                  "'" + SimMeta + "' has correct MetaName",
                  "'" + SimMeta + "' has wrong MetaName");

      printResult(MetaVal == Count, true,
                  "'" + SimMeta + "' has correct MetaVal",
                  "'" + SimMeta + "' has wrong MetaVal");

      // do whatever needed with metadata
      // metaVal can be cast into the appropriate type using
      //   std::any_cast<type>(metaVal)
      Count++;
   }

   // test create() 3
   std::vector<std::shared_ptr<MetaDim>> Dimensions;
   Dimensions.push_back(MetaDim::create(DimName, DimValue));
   const I4 FILL_VALUE = 0;

   auto Data3 = ArrayMetaData::create(
       ArrayName,
       "Description",                   /// long Name or description
       "Units",                         /// units
       "StdName",                       /// CF standard Name
       std::numeric_limits<int>::min(), /// min valid value
       std::numeric_limits<int>::max(), /// max valid value
       FILL_VALUE,                      /// scalar used for undefined entries
       1,                               /// number of dimensions
       Dimensions                       /// dim pointers
   );

   printResult(MetaData::has(ArrayName), true, "'" + ArrayName + "' is created",
               "'" + ArrayName + "' should exist");

   // test get()
   auto Data4 = MetaData::get(ArrayName);

   printResult(Data3 == Data4, true, "get() returns correct instance.",
               "get() returns incorrect instance.");

   // test hasMetaData()
   printResult(Data4->hasMetaData("FillValue"), true,
               "'" + ArrayName + "' has a fill value.",
               "'" + ArrayName + "' does not have a fill value");

   // test getMetaData()
   I4 FillValue;
   ret = Data3->getMetaData("FillValue", FillValue);

   printResult(ret == 0, true, "getMetaData() returns zero",
               "getMetaData() returns non-zero");

   printResult(FILL_VALUE == FillValue, true,
               "'" + ArrayName + "' has a correct fill value",
               "'" + ArrayName + "' has an incorrect fill value");

   // test addMetaData()
   const R8 NewValue = 2.0;
   ret               = Data3->addMetaData("NewMeta", NewValue);

   printResult(ret == 0, true, "addMetaData() returns zero",
               "addMetaData() returns non-zero");

   R8 R8Value;
   ret = Data3->getMetaData("NewMeta", R8Value);

   printResult(NewValue == R8Value, true, "getMetaData() returns correct value",
               "getMetaData() returns incorrect value");

   // test removeMetaData()
   ret = Data3->removeMetaData("NewMeta");

   printResult(ret == 0, true, "removeMetaData() returns zero",
               "removeMetaData() returns non-zero");

   printResult(Data3->hasMetaData("NewMeta"), false,
               "'NewMeta' is removed correctly", "'NewMeta' is not removed");

   // test destroy()
   ret = MetaData::destroy(SimMeta);

   printResult(ret == 0, true, "destroy() returns zero",
               "destroy() returns non-zero");

   printResult(MetaData::has(SimMeta), false,
               "'" + SimMeta + "' is correclty removed",
               "'" + SimMeta + "' is not removed.");
}

void testMetaGroup() {

   const std::string GroupName{"MyGroup"};
   const std::string FieldName{"MyField"};
   const std::string DimName{"MyDim"};
   const I4 DimValue{1};
   int ret;

   // test has()
   printResult(MetaGroup::has(GroupName), false,
               "'" + GroupName + "' is not created",
               "'" + GroupName + "' should not exist");

   // test create()
   auto Group1 = MetaGroup::create(GroupName);

   printResult(MetaGroup::has(GroupName), true,
               "'" + GroupName + "' is created",
               "'" + GroupName + "' should exist");

   // test get()
   auto Group2 = MetaGroup::get(GroupName);

   printResult(Group1 == Group2, true, "get() returns correct instance.",
               "get() returns incorrect instance.");

   // test hasField()
   printResult(Group1->hasField(FieldName), false,
               "'" + FieldName + "' is not in a group",
               "'" + FieldName + "' is in a group");

   // test addField()
   auto Data1 = MetaData::create(FieldName);

   ret = Group1->addField(FieldName);

   printResult(ret == 0, true, "addField() returns zero.",
               "addField() returns non-zero.");

   printResult(Group1->hasField(FieldName), true,
               "'" + FieldName + "' is in a group",
               "'" + FieldName + "' is not in a group");

   std::map<std::string, std::shared_ptr<MetaData>> *Fields;
   Fields = Group1->getAllFields();

   std::string FieldNamePart;
   std::shared_ptr<MetaData> FieldPart;

   // use std::map functionality to loop through all fields
   for (const auto &Pair : *Fields) {
      FieldNamePart = Pair.first;  // retrieve field name part
      FieldPart     = Pair.second; // retrieve field part

      printResult(FieldNamePart == FieldName, true,
                  "Correct FieldName is returned",
                  "Incorrect FieldName is returned");
   }

   // test getField()
   auto Data2 = Group1->getField(FieldName);

   printResult(Data1 == Data2, true, "getField() returns correct instance.",
               "getField() returns incorrect instance.");

   // test removeField()
   ret = Group1->removeField(FieldName);

   printResult(ret == 0, true, "removeField() returns zero.",
               "removeField() returns non-zero.");

   printResult(Group1->hasField(FieldName), false,
               "'" + FieldName + "' is not in a group",
               "'" + FieldName + "' is in a group");

   // test destroy()
   MetaGroup::destroy(GroupName);

   printResult(MetaGroup::has(GroupName), false,
               "'" + GroupName + "' is destroyed correctly",
               "'" + GroupName + "' is not destroyed");
}

int main(int argc, char **argv) {

   int RetVal = 0;

   try {

      testMetaDim();

      testMetaData();

      testMetaGroup();

   } catch (const std::exception &Ex) {
      std::cout << Ex.what() << ": FAIL" << std::endl;
      RetVal -= -1;

   } catch (...) {
      std::cout << "Unknown: FAIL" << std::endl;
      RetVal -= -1;
   }

   return RetVal;
}
