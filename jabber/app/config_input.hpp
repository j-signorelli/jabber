#ifndef JABBER_APP_CONFIG
#define JABBER_APP_CONFIG

#include "params.hpp"

#include <iostream>
#include <format>
#include <sstream>

namespace jabber
{
namespace app
{

/**
 * @defgroup config_group Input Configuration & File Parsing
 * @{
 *
 * @details Presently, Jabber has either ConfigInput for a developer-facing
 * configuration or TOMLConfigInput for runtime-parsing of a TOML-style config
 * file, which is used for all apps in the app suite. See config_template.toml
 * for guidance on preparing a TOML config file.
 */

/// Complete input configuration.
class ConfigInput
{
protected:

   /// Input base flow parameters.
   BaseFlowParams base_flow_;

   /// Input source parameters.
   std::vector<Source::ParamsVariant> sources_;

   /// Input computation parameters.
   CompParams comp_;

   /// Input preCICE parameters.
   std::optional<PreciceParams> precice_;

public:

   /// Get reference to base flow parameters.
   BaseFlowParams& BaseFlow() { return base_flow_; }

   /// Get const reference to base flow parameters.
   const BaseFlowParams& BaseFlow() const { return base_flow_; }

   /// Get reference to source parameters.
   std::vector<Source::ParamsVariant>& Sources() { return sources_; }

   /// Get const reference to source parameters.
   const std::vector<Source::ParamsVariant>& Sources() const { return sources_; }

   /// Get reference to computation parameters.
   CompParams& Comp() { return comp_; }

   /// Get const reference to computation parameters.
   const CompParams& Comp() const { return comp_; }

   /// Get reference to preCICE parameters.
   std::optional<PreciceParams>& Precice() { return precice_; }

   /// Get const reference to preCICE parameters.
   const std::optional<PreciceParams>& Precice() const { return precice_; }

   /// Print the configured base flow parameters.
   void PrintBaseFlowParams(std::ostream &out) const;

   /// Print the configured source parameters.
   void PrintSourceParams(std::ostream &out) const;

   /// Print the configured computation parameters.
   void PrintCompParams(std::ostream &out) const;

   /// Print the configured preCICE parameters.
   void PrintPreciceParams(std::ostream &out) const;
};

/**
 * @brief TOML config file input.
 * 
 * @details This class accepts a TOML config file and parses it appropriately.
 * Separate static parsing functions are defined which accept a TOML-formatted
 * string for flexibility and simpler unit testing.
 */
class TOMLConfigInput : public ConfigInput
{
public:

   /**
    * @brief Parse BaseFlowParams from a serialized TOML string of that
    * section.
    */
   static void ParseBaseFlow(std::string toml_string, 
                              BaseFlowParams &op);

   /**
    * @brief Parse InputXY parameters from a serialized TOML string of that
    * section.
    */
   static void ParseInputXY(std::string toml_string,
                              InputXY::ParamsVariant &opv);

   /**
    * @brief Parse FunctionType parameters from a serialized TOML string of that
    * section.
    */
   static void ParseFunctionType(std::string toml_string,
                              FunctionType::ParamsVariant &opv);

   /**
    * @brief Parse DiscMethod parameters from a serialized TOML string of
    * that section.
    */
   static void ParseDiscMethod(std::string toml_string,
                                 DiscMethod::ParamsVariant &opv);
   
   /**
    * @brief Parse Direction parameters from a serialized TOML string of that
    * section.
    */
   static void ParseDirection(std::string toml_string,
                              Direction::ParamsVariant &opv);

   /**
    * @brief Parse TransferFunction parameters from a serialized TOML string of
    * that section.
    */
   static void ParseTransferFunction(std::string toml_string,
                              TransferFunction::ParamsVariant &opv);

   /**
    * @brief Parse Source parameters from a serialized TOML string of that
    * section.
    */
   static void ParseSource(std::string toml_string,
                           Source::ParamsVariant &opv);

   /**
    * @brief Parse Computation parameters from a serialized TOML string of that
    * section.
    */
   static void ParseComputation(std::string toml_string,
                                 CompParams &op);

   /**
    * @brief Parse PreciceParams from a serialized TOML string of that
    * section.
    */
   static void ParsePrecice(std::string toml_string,
                              PreciceParams &op);


   /// Construct an uninitialized TOMLConfigInput object.
   TOMLConfigInput() {};

   /**
    * @brief Construct a new TOMLConfigInput object.
    * 
    * @param config_file      TOML config file address to parse.
    * @param out              [Optional] ostream to write parsed config file
    *                         to (verbose processing).
    */
   TOMLConfigInput(std::string config_file, std::ostream *out=nullptr);

};

/// @}
// end of config_group

} // namespace app

} // namespace jabber

#endif // JABBER_APP_CONFIG
