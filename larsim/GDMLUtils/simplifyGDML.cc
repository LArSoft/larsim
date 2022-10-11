/**
 * @file   simplifyGDML.cc
 * @brief  Reprocesses a GDML file via GEANT4.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 13, 2017
 *
 * Run with `--help` argument for usage instructions.
 *
 * This program needs to be linked to:
 * - GEANT4
 * - CLHEP
 * - XERCES (C)
 *
 * Example of build command in UPS environment:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * g++ -Wall -pedantic -std=c++14 \
 *   -I"${GEANT4_FQ_DIR}/include" -I"$CLHEP_INC" -I"$XERCES_C_INC" \
 *   -L"${GEANT4_FQ_DIR}/lib64" -lG4persistency -lG4geometry \
 *   -L"$XERCES_C_LIB" -lxerces-c \
 *   -o simplifyGDML.exe simplifyGDML.cc
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Note that this program is a glorified 3-line code: parse XML file, get
 * the world volume, write it.
 *
 */

// GEANT4
#include "Geant4/G4GDMLParser.hh"

// POSIX/UNIX
#include <getopt.h> // getopt_long(), option
#include <unistd.h> // access()

// C/C++ standard libraries
#include <cstdio>  // std::remove()
#include <cstdlib> // std::exit()
#include <iostream>
#include <string>

//------------------------------------------------------------------------------
// argument parsing
struct ConfigurationParameters {

  static const unsigned int StartingDebugLevel = 0;
  static std::string const DefaultSetupName;
  static std::string const DefaultSchemaURL;

  std::string sourcePath; ///< GDML input file name.
  std::string destPath;   ///< GDML output file name.

  /// Name of the chosen setup in the GDML source.
  std::string setupName = DefaultSetupName;

  /// URL of the schema used while writing GDML.
  std::string schemaPath = DefaultSchemaURL;

  bool validate = false;                        ///< Ask Geant4 to validate the source.
  bool overwrite = false;                       ///< Overwrite the output file if already present.
  bool dontWrite = false;                       ///< Do not write the output file.
  bool help = false;                            ///< Print usage instructions and exit.
  unsigned int debugLevel = StartingDebugLevel; ///< Debug level.

  void print() const;

}; // struct ConfigurationParameters

std::string const ConfigurationParameters::DefaultSetupName = "Default";
std::string const ConfigurationParameters::DefaultSchemaURL = G4GDML_DEFAULT_SCHEMALOCATION;

void ConfigurationParameters::print() const
{

  std::cout << "Configuration:"
               "\n  input file:  '"
            << sourcePath
            << "'"
               "\n  output file: '"
            << destPath
            << "'"
               "\n  setup name:  '"
            << setupName
            << "'"
               "\n  schema path: '"
            << schemaPath
            << "'"
               "\n  validate:     "
            << std::boolalpha << validate << "\n  only read:    " << std::boolalpha << dontWrite
            << "\n  overwrite:    " << std::boolalpha << overwrite
            << "\n  print help:   " << std::boolalpha << help << "\n  debug level:  " << debugLevel
            << std::endl;

} // ConfigurationParameters::print()

struct ConfigurationParser {

  static const unsigned int DefaultDebugLevel = 1;

  typedef enum {
    Success,
    InvalidNumber,
    InvalidOption,
    LogicError,
    ExtraArguments,
    MissingArgument
  } error;

  static error parse(ConfigurationParameters& params, unsigned int argc, char** argv);

  static void printHelp(int exitCode = 0, const char* progName = "simplifyGDML");

private:
  struct opt; // single letter options namespace

}; // struct ConfigurationParameters

void ConfigurationParser::printHelp(int exitCode /* = 0 */,
                                    const char* progName /* = "simplifyGDML" */)
{
  std::cout << "Asks GEANT4 to process and rewrite a GDML file."
               " This effectively simplifies the complexity of constructs in the GDML"
               " code."
               "\n"
               "\nUsage:  "
            << progName
            << "  [options] [--] inputPath [outputPath]"
               "\n"
               "\nThe output file will contain the world volume of the specified setup."
               "\nIt must be different from the input file:"
               " GEANT4 will refuse to overwrite."
               "\n"
               "\nNOTE: the path to the GDML schema in the output may need to be fixed by"
               " hand."
               "\nTo allow validation, the GDML schema must be present as described in the"
               " header of the GDML file."
               "\nLArSoft does not necessarily distributes GDML schema, so these lines"
               " are aimed to load them from network:"
               "\n"
               "\n<gdml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
               "\n      "
               "xsi:noNamespaceSchemaLocation=\"http://service-spi.web.cern.ch/service-spi/app/"
               "releases/GDML/schema/gdml.xsd\""
               "\n      >"
               "\n"
               "\nOptions:"
               "\n--validate , -V"
               "\n    ask Geant4 to validated the GDML while reading (validation output"
               "\n    will be on screen, with no effect to the rest of the program)"
               "\n--overwrite , -f"
               "\n    not implemented yet"
               "\n--nowrite , -r"
               "\n    only read (and validate as requested), do not write a new GDML file"
               "\n--setup=SETUPNAME , -s SETUPNAME"
               "\n    name of path or URL to the schema to use for writing."
               "\n    By default, the '"
            << ConfigurationParameters::DefaultSetupName
            << "' setup is chosen."
               "\n--schema=SCHEMAURL , -S SCHEMAURL"
               "\n    path or URL to the schema to use for writing."
               "\n    By default, the '"
            << ConfigurationParameters::DefaultSchemaURL
            << "' schema is used."
               "\n--help , -h , -?"
               "\n    print these usage instructions and exit"
               "\n"
            << std::endl;
  if (exitCode >= 0) std::exit(exitCode);
} // ConfigurationParser::printHelp()

struct ConfigurationParser::opt {
  static constexpr const char Validate = 'v';
  static constexpr const char Schema = 'S';
  static constexpr const char Setup = 's';
  static constexpr const char NoWrite = 'r';
  static constexpr const char Overwrite = 'f';
  static constexpr const char Debug = 'd';
  static constexpr const char Help = 'h';
}; // ConfigurationParser::opt

ConfigurationParser::error ConfigurationParser::parse(ConfigurationParameters& params,
                                                      unsigned int argc,
                                                      char** argv)
{
  static const option longopts[] = {{"validate", no_argument, NULL, opt::Validate},
                                    {"schema", required_argument, NULL, opt::Schema},
                                    {"setup", required_argument, NULL, opt::Setup},
                                    {"overwrite", no_argument, NULL, opt::Overwrite},
                                    {"nowrite", no_argument, NULL, opt::NoWrite},
                                    {"debug", optional_argument, NULL, opt::Debug},
                                    {"help", no_argument, NULL, opt::Help},
                                    //    { "one",      no_argument,       &value, 1 },
                                    //    { "two",      no_argument,       &value, 2 },
                                    {NULL, 0, NULL, 0}}; // longopts

  // automatically build the short option string from the long one
  std::string shortopts = ":"; // this means no error printout
  for (auto const& longopt : longopts) {
    if (!longopt.name) break;
    if (longopt.val == 0) continue;
    shortopts += (char)longopt.val;
    if (longopt.has_arg != no_argument) shortopts += ':';
    if (longopt.has_arg == optional_argument) shortopts += ':';
  } // for

  // ----------------------------------------------------------------------
  // options
  error res = error::Success;
  char ch;
  optind = 1;
  bool forgiving = false;
  while ((ch = getopt_long(argc, argv, shortopts.c_str(), longopts, NULL)) != -1) {
    switch (ch) {
    case opt::Validate: params.validate = true; continue;
    case opt::Setup: params.setupName = optarg; continue;
    case opt::Schema: params.schemaPath = optarg; continue;
    case opt::Overwrite: params.overwrite = true; continue;
    case opt::NoWrite: params.dontWrite = true; continue;
    case opt::Debug:
      if (optarg) {
        std::istringstream sstr(optarg);
        sstr >> params.debugLevel;
        if (!sstr) {
          std::cerr << "Invalid debug level: '" << optarg << "'" << std::endl;
          res = error::InvalidNumber;
          continue;
        }
      }
      else
        params.debugLevel = DefaultDebugLevel;
      continue;
    case '?': // this is special since has a specific meaning for getopt()
      if (optopt != '?') {
        if (!forgiving) {
          std::cerr << "Invalid option: '" << argv[optind - 1] << "'" << std::endl;
        }
        res = error::InvalidOption;
        continue;
      }
      // follow into:
    case opt::Help:
      params.help = true;
      forgiving = true;
      continue;
    default:
      std::cerr << "Internal error: option '" << ch << "' not supported yet." << std::endl;
      res = error::LogicError;
      continue;
    } // switch
  }   // while

  // ----------------------------------------------------------------------
  // arguments
  if (optind >= (int)argc) {
    if (!forgiving) std::cerr << "Source file name is required!" << std::endl;
    return error::MissingArgument;
  }
  params.sourcePath = argv[optind++];
  if (optind < (int)argc) params.destPath = argv[optind++];

  if (optind < (int)argc) {
    if (!forgiving) {
      std::cerr << "Spurious arguments: '" << argv[optind] << "'";
      if (optind + 1 < (int)argc) std::cerr << " and " << (argc - optind - 1) << " more";
      std::cerr << "." << std::endl;
    }
    return error::ExtraArguments;
  }

  // ----------------------------------------------------------------------
  return res;
} // ConfigurationParser::parse()

//------------------------------------------------------------------------------
std::string addNameSuffix(std::string name, std::string suffix)
{
  auto const iExt = name.find(".gdml");
  if (iExt == std::string::npos)
    name.append(suffix);
  else
    name.insert(iExt, suffix);
  return name;
} // addNameSuffix()

bool exists(std::string path)
{
  return access(path.c_str(), F_OK);
} // exists()

//------------------------------------------------------------------------------
G4VPhysicalVolume* readWorldVolume(G4GDMLParser& parser, ConfigurationParameters const& params)
{
  decltype(auto) worldVolume = parser.GetWorldVolume(params.setupName);
  if (!worldVolume) {
    std::cerr << "Failed to find the world volume for setup '" << params.setupName << "'"
              << std::endl;
  }
  return worldVolume;
} // readWorldVolume()

int writeWorld(G4GDMLParser& parser, ConfigurationParameters const& params)
{

  std::cout << std::string(80, '-') << "\nFetching the world volume for setup '" << params.setupName
            << "'..." << std::endl;
  decltype(auto) worldVolume = readWorldVolume(parser, params);
  if (!worldVolume) {
    std::cerr << "Failed to find the world volume for setup '" << params.setupName << "'"
              << std::endl;
    return 1;
  }

  if (exists(params.destPath)) {
    if (params.overwrite) {
      std::cout << "Overwriting the destination file..." << std::endl;
      std::remove(params.destPath.c_str());
    }
    else {
      std::cerr << "Destination file '" << params.destPath << "' already exists." << std::endl;
      return 2;
    }
  }

  std::cout << std::string(80, '-') << "\nWriting '" << params.destPath << "'..." << std::endl;
  parser.Write(params.destPath, worldVolume, false /* do not add references to names*/);

  std::cout << std::string(80, '-') << "\n'" << params.sourcePath << "' => [" << params.setupName
            << "] => '" << params.destPath << "': done." << std::endl;

  return 0;
} // writeVolume()

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{

  //
  // argument parsing
  //
  ConfigurationParameters params;

  auto parseRes = ConfigurationParser::parse(params, argc, argv);

  if (params.destPath.empty() && !params.dontWrite)
    params.destPath = addNameSuffix(params.sourcePath, "-simplified");

  if (params.debugLevel > 0) params.print();

  // print usage instructions before exiting
  if (params.help) {
    ConfigurationParser::printHelp(0, argv[0]);
    return 0;
  }
  if (parseRes != ConfigurationParser::Success) return (int)parseRes;

  //
  // the magic
  //
  G4GDMLParser parser;

  std::cout << std::string(80, '-') << "\nReading '" << params.sourcePath << "'..." << std::endl;
  parser.Read(params.sourcePath, params.validate);

  std::cout << std::string(80, '-') << "\nFetching the world volume for setup '" << params.setupName
            << "'..." << std::endl;
  // for c2: worldVolume is unused here
  //decltype(auto) worldVolume = parser.GetWorldVolume(params.setupName);
  parser.GetWorldVolume(params.setupName);

  if (params.dontWrite) {
    if (!params.destPath.empty()) {
      std::cerr << "Output path ignored since we don't write anything." << std::endl;
    }
  }
  else {
    int res = writeWorld(parser, params);
    if (res != 0) return res;
  }

  //
  // the missing error code
  //
  return 0;
} // main()
